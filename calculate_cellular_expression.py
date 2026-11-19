import argparse
import pandas as pd
import os
import time
import numpy as np
from scipy import ndimage
from glob import glob
import pickle
from tqdm import tqdm
from treelib import Tree
import multiprocessing as mp

from Utils.data_io import nib_load, check_folder
from Utils.data_structure import read_cd_file
from lineage_stat.data_structure import construct_celltree



def save_expression_to_readable_file(name_dictionary_path, embryo_name, running_max_time, CD_file_root, stat_file_root):
    label_name_dict = pd.read_csv(name_dictionary_path, index_col=0).to_dict()['0']
    # name_label_dict = {value: key for key, value in label_name_dict.items()}

    # ---------------------read cd file---------------------------------
    cd_file_path = os.path.join(CD_file_root, 'CD{}.csv'.format(embryo_name))
    # cd_file_all_dict = {}
    cd_file_dataframe = read_cd_file(cd_file_path)
    # print('constructing cell tree for ', embryo_name)
    cell_tree = construct_celltree(cd_file_path, running_max_time, name_dictionary_path)
    # ---------------------------------------------------------------

    # -==============assemble ========== surface-------------and-----------volume===================================
    all_names = [cname for cname in cell_tree.expand_tree(mode=Tree.WIDTH)]
    # for idx, cell_name in enumerate(all_names):
    volume_embryo = pd.DataFrame(
        np.full(shape=(running_max_time, len(all_names)), fill_value=np.nan, dtype=np.float32),
        index=range(1, running_max_time + 1), columns=all_names)

    for tp in tqdm(range(1, running_max_time + 1),
                   desc='assembling volume and surface area of {} result'.format(embryo_name)):
        path_tmp = os.path.join(stat_file_root, embryo_name)
        with open(os.path.join(path_tmp, '{}_{}_cellwise_expression.txt'.format(embryo_name, str(tp).zfill(3))),
                  'rb') as handle:
            cellwise_expression_dict = pickle.load(handle)

        for cell_label_, vol_value in cellwise_expression_dict.items():
            cell_name_ = label_name_dict[cell_label_]
            volume_embryo.loc[tp, cell_name_] = vol_value.astype(np.float32)

    volume_embryo = volume_embryo.loc[:, ((volume_embryo != 0) & (~np.isnan(volume_embryo))).any(axis=0)]
    volume_embryo.to_csv(os.path.join(stat_file_root, embryo_name + '_cellwise_expression.csv'))

    # ---------------------assemble surface area expression---------------------------------
    surface_embryo = pd.DataFrame(
        np.full(shape=(running_max_time, len(all_names)), fill_value=np.nan, dtype=np.float32),
        index=range(1, running_max_time + 1), columns=all_names)

    for tp in tqdm(range(1, running_max_time + 1),
                   desc='assembling surface area of {} result'.format(embryo_name)):
        path_tmp = os.path.join(stat_file_root, embryo_name)
        with open(os.path.join(path_tmp, '{}_{}_surfacewise_expression.txt'.format(embryo_name, str(tp).zfill(3))),
                  'rb') as handle:
            surfacewise_expression_dict = pickle.load(handle)

        for cell_label_, surface_value in surfacewise_expression_dict.items():
            cell_name_ = label_name_dict[cell_label_]
            surface_embryo.loc[tp, cell_name_] = surface_value.astype(np.float32)

    surface_embryo = surface_embryo.loc[:, ((surface_embryo != 0) & (~np.isnan(surface_embryo))).any(axis=0)]
    surface_embryo.to_csv(os.path.join(stat_file_root, embryo_name + '_surfacewise_expression.csv'))

    # ------------contact-----------initialize the contact csv  file----------------------
    # Get tuble lists with elements from the list

    name_combination = []
    first_level_names = []
    for i, name1 in enumerate(all_names):
        for name2 in all_names[i + 1:]:
            if not (cell_tree.is_ancestor(name1, name2) or cell_tree.is_ancestor(name2, name1)):
                first_level_names.append(name1)
                name_combination.append((name1, name2))

    multi_index = pd.MultiIndex.from_tuples(name_combination, names=['cell1', 'cell2'])
    # print(multi_index)
    stat_embryo = pd.DataFrame(
        np.full(shape=(running_max_time, len(name_combination)), fill_value=np.nan, dtype=np.float32),
        index=range(1, running_max_time + 1), columns=multi_index)
    # set zero element to express the exist of the specific nucleus
    for cell_name in all_names:
        if cell_name not in first_level_names:
            continue
        try:
            cell_time = cell_tree.get_node(cell_name).data.get_time()
            cell_time = [x for x in cell_time if x <= running_max_time]
            stat_embryo.loc[cell_time, (cell_name, slice(None))] = 0
        except:
            cell_name
    # print(stat_embryo)
    # edges_view = point_embryo.edges(data=True)
    for tp in tqdm(range(1, running_max_time + 1),
                   desc='assembling contact surface of {} result'.format(embryo_name)):
        path_tmp = os.path.join(stat_file_root, embryo_name)
        with open(os.path.join(path_tmp, '{}_{}_contactwise_expression.txt'.format(embryo_name, str(tp).zfill(3))),
                  'rb') as handle:
            contact_expression_dict = pickle.load(handle)
        for contact_sur_idx, contact_sur_value in contact_expression_dict.items():
            [cell1, cell2] = contact_sur_idx.split('_')
            cell1_name = label_name_dict[int(cell1)]
            cell2_name = label_name_dict[int(cell2)]
            if (cell1_name, cell2_name) in stat_embryo.columns:
                stat_embryo.loc[tp, (cell1_name, cell2_name)] = contact_sur_value.astype(np.float32)
            elif (cell2_name, cell1_name) in stat_embryo.columns:
                stat_embryo.loc[tp, (cell2_name, cell1_name)] = contact_sur_value.astype(np.float32)
            else:
                pass
                # print('columns missing (cell1_name, cell2_name)')
    stat_embryo = stat_embryo.loc[:, ((stat_embryo != 0) & (~np.isnan(stat_embryo))).any(axis=0)]
    # print(stat_embryo)
    stat_embryo.to_csv(os.path.join(stat_file_root, embryo_name + '_contactwise_expression.csv'))
    # --------------------------------------------------------------------------------------------


def calculate_volume_surface_contactwise_expression(config):
    embryo_seg_path = config[0]
    embryo_expression_root_path = config[1]
    saving_expression_data = config[2]

    embryo_name, time_point = os.path.basename(embryo_seg_path).split('_')[:2]
    embryo_expression_path = os.path.join(embryo_expression_root_path, embryo_name, 'Gene',
                                          '{}_{}_expression.nii.gz'.format(embryo_name, time_point))

    if embryo_seg_path is None:
        raise Exception("Sorry, no path_embryo given!")

    # ------------------------calculate surface points using dialation for each cell --------------------
    frame_this_embryo = str(time_point).zfill(3)
    file_name = embryo_name + '_' + frame_this_embryo + '.nii.gz'

    segmented_arr = nib_load(embryo_seg_path).astype(int)
    raw_expression_arr = nib_load(embryo_expression_path)
    # .transpose([2, 1, 0])

    start_time = time.time()

    # ======================cell contact pairs==========================

    # -----!!!!!!!!get cell edge boundary!!!!--------------
    cell_bone_mask = np.zeros(segmented_arr.shape)
    # start_time=time.time()
    for cell_idx in np.unique(segmented_arr)[1:]:
        this_cell_arr = (segmented_arr == cell_idx)
        this_cell_arr_dilation = ndimage.binary_dilation(this_cell_arr, iterations=1)
        cell_bone_mask[np.logical_xor(this_cell_arr_dilation, this_cell_arr)] = 1
    print('edge boundary timing', time.time() - start_time)

    # cell_mask = segmented_arr != 0
    # boundary_mask_tmp = np.logical_xor(cell_mask == 0 , ndimage.binary_dilation(cell_mask))
    # boundary_mask=np.zeros(segmented_arr.shape)
    # boundary_mask[boundary_mask_tmp]=1
    #
    # nib_save(cell_bone_mask,r'./test_1.nii.gz')
    # print((cell_bone_mask==1).sum())
    # nib_save(boundary_mask,r'./test_2.nii.gz')
    # print((boundary_mask==1).sum())

    [x_bound, y_bound, z_bound] = np.nonzero(cell_bone_mask)
    boundary_elements = []

    # find boundary/points between cells
    start_time = time.time()
    for (x, y, z) in zip(x_bound, y_bound, z_bound):
        neighbors = segmented_arr[np.ix_(range(x - 1, x + 2), range(y - 1, y + 2), range(z - 1, z + 2))]
        neighbor_labels = list(np.unique(neighbors))[1:]  # with order
        # print(neighbor_labels)
        # if 0 in neighbor_labels:
        #     neighbor_labels.remove(0)
        if len(neighbor_labels) == 2:  # contact between two cells
            boundary_elements.append(sorted(neighbor_labels))
    print('getting contact pairs', time.time() - start_time)

    cell_contact_pairs = list(np.unique(np.array(boundary_elements), axis=0))
    # ===================================================================

    volume_dict = {}
    surface_dict = {}
    contact_expression_dict = {}

    # =================contact wise==============================================
    start_time = time.time()
    for (label1, label2) in cell_contact_pairs:
        contact_mask_tmp = np.logical_and(ndimage.binary_dilation(segmented_arr == label1),
                                          ndimage.binary_dilation(segmented_arr == label2))
        contact_mask = np.logical_and(contact_mask_tmp, cell_bone_mask)
        contact_sum = contact_mask.sum()
        dilated_contact_mask = ndimage.binary_dilation(contact_mask)
        # print(dilated_contact_mask.sum())
        expression_sum = sum(raw_expression_arr[dilated_contact_mask])
        if contact_sum > 2:
            # cell_contact_pair_renew.append((label1, label2))
            str_key = str(label1) + '_' + str(label2)
            # contact_area_dict[str_key] = 0
            contact_expression_dict[str_key] = expression_sum
            print(str_key, expression_sum)
    print('get pair contact surface expression time', time.time() - start_time)
    # =============================================================================

    cell_list = np.unique(segmented_arr)

    for cell_key in cell_list:
        if cell_key != 0:
            dilated_cell_mask = ndimage.binary_dilation(segmented_arr == cell_key)
            expression_this_cell = sum(raw_expression_arr[dilated_cell_mask])
            volume_dict[cell_key] = expression_this_cell
            print(cell_key, expression_this_cell)

            this_cell_surface_mask = ndimage.binary_dilation(
                np.logical_xor(dilated_cell_mask, (segmented_arr == cell_key)))
            expression_this_cell_surface = sum(raw_expression_arr[this_cell_surface_mask])
            surface_dict[cell_key] = expression_this_cell_surface

    path_tmp = os.path.join(saving_expression_data, embryo_name)
    check_folder(path_tmp)
    with open(os.path.join(path_tmp, file_name.split('.')[0] + '_cellwise_expression.txt'), 'wb+') as handle:
        pickle.dump(volume_dict, handle, protocol=4)
    with open(os.path.join(path_tmp, file_name.split('.')[0] + '_surfacewise_expression.txt'), 'wb+') as handle:
        pickle.dump(surface_dict, handle, protocol=4)
    with open(os.path.join(path_tmp, file_name.split('.')[0] + '_contactwise_expression.txt'), 'wb+') as handle:
        pickle.dump(contact_expression_dict, handle, protocol=4)

    # -------------------------------------------------------------------------------------------------------


def save_expression_by_label(embryo_name, running_max_time, stat_file_root):
    """Assemble the per-time-point results into tables whose columns are segmentation labels.

    Used when no lineage file / name dictionary is given, so cell names are not required.
    """
    path_tmp = os.path.join(stat_file_root, embryo_name)
    for kind in ('cellwise', 'surfacewise', 'contactwise'):
        rows = {}
        for tp in range(1, running_max_time + 1):
            file_path = os.path.join(path_tmp, '{}_{}_{}_expression.txt'.format(embryo_name, str(tp).zfill(3), kind))
            with open(file_path, 'rb') as handle:
                rows[tp] = {str(key): float(value) for key, value in pickle.load(handle).items()}
        table = pd.DataFrame.from_dict(rows, orient='index')
        table = table[sorted(table.columns, key=lambda c: [int(x) for x in c.split('_')])]
        table.index.name = 'time_point'
        table.to_csv(os.path.join(stat_file_root, '{}_{}_expression_by_label.csv'.format(embryo_name, kind)))


def parse_args():
    parser = argparse.ArgumentParser(
        description='Per-cell, cell-surface and cell-contact fluorescence quantification.')
    parser.add_argument('--seg-root', required=True,
                        help='Folder with <embryo>/<seg-dir-name>/<embryo>_<tp>_segCell.nii.gz')
    parser.add_argument('--seg-dir-name', default='SegCell',
                        help='Per-embryo folder holding the named-cell segmentations (default SegCell)')
    parser.add_argument('--expression-root', required=True,
                        help='Folder with <embryo>/Gene/<embryo>_<tp>_expression.nii.gz')
    parser.add_argument('--cd-root', default=None,
                        help='Optional. Folder with CD<embryo>.csv (lineage file), for cell-name columns')
    parser.add_argument('--name-dictionary', default=None,
                        help='Optional. name_dictionary.csv mapping segmentation labels to cell names')
    parser.add_argument('--output', required=True, help='Output folder for the expression tables')
    parser.add_argument('--embryo', action='append',
                        help='Embryo name; may be repeated. Default: every embryo under --seg-root')
    parser.add_argument('--workers', type=int, default=1, help='Number of parallel processes (default 1)')
    parser.add_argument('--assemble-only', action='store_true',
                        help='Skip the calculation; only assemble existing intermediate files into CSV tables')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    if (args.cd_root is None) != (args.name_dictionary is None):
        raise SystemExit('Give both --cd-root and --name-dictionary (cell-name tables), or neither (label tables).')
    use_cell_names = args.cd_root is not None
    if use_cell_names and not os.path.isfile(args.name_dictionary):
        raise SystemExit('Name dictionary not found: {}'.format(args.name_dictionary))
    os.makedirs(args.output, exist_ok=True)

    embryo_names = args.embryo or sorted(
        d for d in os.listdir(args.seg_root)
        if os.path.isdir(os.path.join(args.seg_root, d, args.seg_dir_name)))
    if not embryo_names:
        raise SystemExit('No <embryo>/{} folders found under {}'.format(args.seg_dir_name, args.seg_root))

    for embryo_name in embryo_names:
        seg_files = sorted(glob(os.path.join(args.seg_root, embryo_name, args.seg_dir_name, '*.nii.gz')))
        if not seg_files:
            print('skip {} (no segmentation files)'.format(embryo_name))
            continue
        max_time = len(seg_files)
        if use_cell_names:
            cd_file = os.path.join(args.cd_root, 'CD{}.csv'.format(embryo_name))
            if not os.path.isfile(cd_file):
                raise SystemExit('Lineage file not found: {}'.format(cd_file))

        if not args.assemble_only:
            # check every expression image exists before the long calculation starts
            missing = []
            for seg_file in seg_files:
                tp = os.path.basename(seg_file).split('_')[1]
                expression_file = os.path.join(args.expression_root, embryo_name, 'Gene',
                                               '{}_{}_expression.nii.gz'.format(embryo_name, tp))
                if not os.path.isfile(expression_file):
                    missing.append(expression_file)
            if missing:
                raise SystemExit('{} expression image(s) missing, e.g. {}'.format(len(missing), missing[0]))

            configs = [[seg_file, args.expression_root, args.output] for seg_file in seg_files]
            desc = 'calculating {} expression'.format(embryo_name)
            if args.workers > 1:
                with mp.Pool(min(args.workers, len(configs))) as pool:
                    for _ in tqdm(pool.imap_unordered(calculate_volume_surface_contactwise_expression, configs),
                                  total=len(configs), desc=desc):
                        pass
            else:
                for config in tqdm(configs, desc=desc):
                    calculate_volume_surface_contactwise_expression(config)

        # turn the per-time-point results into three tables (rows: time points)
        if use_cell_names:
            save_expression_to_readable_file(args.name_dictionary, embryo_name, max_time, args.cd_root, args.output)
        else:
            save_expression_by_label(embryo_name, max_time, args.output)
        print('{}: tables written to {}'.format(embryo_name, args.output))
