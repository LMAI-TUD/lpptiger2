# -*- coding: utf-8 -*-
#
# Copyright (C) 2021-2022  LMAI_team @ TU Dresden:
#     LMAI_team: Zhixu Ni, Maria Fedorova
# Copyright (C) 2016-2021  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
#     SysMedOs_team: Zhixu Ni, Georgia Angelidou, Maria Fedorova
#
# LPPtiger2 is Dual-licensed
#     For academic and non-commercial use: `AGPL License V3` Please read more information by the following link:
#         [GNU Affero General Public License] (https://www.gnu.org/licenses/agpl-3.0.en.html)
#     For commercial use:
#         Please contact the LMAI_team by email.
#         https://tu-dresden.de/med/mf/zml/forschungsgruppen/fedorova/
#
# Please cite our publication in an appropriate form:
#     Ni, Zhixu, Georgia Angelidou, Ralf Hoffmann, and Maria Fedorova.
#     LPPtiger software for lipidome-specific prediction
#     and hunter of oxidized phospholipids from LC-MS datasets
#     Scientific Reports 7, Article number: 15138 (2017).
#     DOI: 10.1038/s41598-017-15363-z
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@tu-dresden.de
#     LPPtiger2 repository: https://github.com/LMAI-TUD/lpptiger2
#
import re

from natsort import natsorted
import pandas as pd

from libtiger.utils.default_values import charge_cfg_dct


# def filter_res_adducts(lipid_dct: dict, pr_adduct: str):
#
#     pr_dct = {}
#
#     pr_info_dct = lipid_dct.copy()
#
#     pr_name = lipid_dct.get("name", "")
#     lipid_info_dct = lipid_dct.get("info", {}).copy()
#     lipid_class = lipid_info_dct.get("lipid_class", "")
#
#     is_valid_pr_adduct = False
#
#     if lipid_class.endswith("PC"):
#         if re.findall(r"COOH", pr_name):
#             if pr_adduct == "[M-H]-":
#                 is_valid_pr_adduct = True
#             else:
#                 pass
#         else:
#             if pr_adduct in ["[M+HCOO]-", "[M+CH3COO]-"]:
#                 is_valid_pr_adduct = True
#             else:
#                 pass
#     else:
#         is_valid_pr_adduct = True
#     if is_valid_pr_adduct:
#         res_lst = []
#         for fa in ["FA1", "FA2", "FA3", "FA4"]:
#             if fa in lipid_info_dct:
#                 res_lst.append(fa)
#         lipid_class_chg_adduct = charge_cfg_dct.get(lipid_class, {}).get(pr_adduct, [])
#         for res in res_lst:
#             res_mz_info = lipid_info_dct[res].get("mz_info", {}).copy()
#             filtered_res_mz_info = {}
#             for res_adduct in res_mz_info:
#                 if res_adduct in lipid_class_chg_adduct:
#                     filtered_res_mz_info[res_adduct] = res_mz_info[res_adduct]
#                 else:
#                     pass
#             lipid_info_dct[res]["mz_info"] = filtered_res_mz_info
#         pr_info_dct["info"] = lipid_info_dct
#
#         pr_dct = {f"{pr_name}#{pr_adduct}": pr_info_dct}
#
#     return pr_dct


def build_search_space(epilipidome: dict, charge_modes: list, adducts: list, mz_range: list = None):

    raw_lipid_space_dct = {}
    raw_residue_class_dct = {}
    raw_residue_frag_dct = {}
    raw_residue_nl_dct = {}

    if mz_range is None:
        mz_min = 200
        mz_max = 1200
    elif len(mz_range) == 2:
        mz_min = mz_range[0]
        mz_max = mz_range[1]
    else:
        mz_min = 200
        mz_max = 1200

    space_adducts_info = {}

    for charge_mode in charge_modes:
        if re.match(r"^neg.*$", charge_mode, re.IGNORECASE):
            neg_pr_adducts = ["[M-H]-"]
            if "[M+HCOO]-" in adducts:
                neg_pr_adducts.append("[M+HCOO]-")
            if "[M+CH3COO]-" in adducts:
                neg_pr_adducts.append("[M+CH3COO]-")
            if "[M+HCOO]-" not in adducts and "[M+CH3COO]-" not in adducts:
                neg_pr_adducts.append("[M+HCOO]-")
            space_adducts_info["neg"] = list(set(neg_pr_adducts))
        elif re.match(r"^pos.*$", charge_mode, re.IGNORECASE):
            pos_pr_adducts = []
            for pos_adduct in ["[M+H]+", "[M+NH4]+", "[M+Na]+"]:
                if pos_adduct in adducts:
                    pos_pr_adducts.append(pos_adduct)
            if not pos_pr_adducts:
                pos_pr_adducts = ["[M+H]+", "[M+NH4]+", "[M+Na]+"]
            space_adducts_info["pos"] = list(set(pos_pr_adducts))
        else:
            pass

    if not space_adducts_info:
        space_adducts_info = {
            "neg": ["[M-H]-", "[M+HCOO]-", "[M+CH3COO]-"],
            "pos": ["[M+H]+", "[M+NH4]+", "[M+Na]+"]
        }

    adduct_lst = space_adducts_info.get("neg", []) + space_adducts_info.get("pos", [])

    for lipid_name in epilipidome:
        lipid_adduct_dct = epilipidome[lipid_name]
        for pr_label in lipid_adduct_dct:
            lipid_dct = lipid_adduct_dct[pr_label]
            pr_adduct = lipid_dct.get("adduct")
            pr_mz = float(lipid_dct.get("mz", 0.0))
            if pr_adduct in adduct_lst and mz_min <= pr_mz <= mz_max:
                chg_formula = lipid_dct.get("formula", "")
                if chg_formula not in raw_lipid_space_dct:
                    raw_lipid_space_dct[chg_formula] = {pr_label: lipid_dct}
                else:
                    raw_lipid_space_dct[chg_formula].update({pr_label: lipid_dct})

                ions_dct = lipid_dct.get("ions", {})
                class_dct = ions_dct.get("class_info", {})
                frag_dct = ions_dct.get("frag_info", {})
                nl_dct = ions_dct.get("nl_info", {})
                
                for class_label in class_dct:
                    class_formula = class_dct[class_label].get("formula", "")
                    if class_formula not in raw_residue_class_dct:
                        raw_residue_class_dct[class_formula] = {class_label: class_dct[class_label]}
                    else:
                        raw_residue_class_dct[class_formula].update({class_label: class_dct[class_label]})

                for frag_label in frag_dct:
                    # frag_mz = frag_dct[frag_label].get("mz", 0.0)
                    # if frag_mz not in raw_residue_frag_dct:
                    #     raw_residue_frag_dct[frag_mz] = [frag_label]
                    # else:
                    #     raw_residue_frag_dct[frag_mz].append(frag_label)
                    frag_formula = frag_dct[frag_label].get("formula", "")
                    if frag_formula not in raw_residue_frag_dct:
                        raw_residue_frag_dct[frag_formula] = {frag_label: frag_dct[frag_label]}
                    else:
                        raw_residue_frag_dct[frag_formula].update({frag_label: frag_dct[frag_label]})
                for nl_label in nl_dct:
                    # nl_mz = nl_dct[nl_label].get("mz", 0.0)
                    # if nl_mz not in raw_residue_nl_dct:
                    #     raw_residue_nl_dct[nl_mz] = [nl_label]
                    # else:
                    #     raw_residue_nl_dct[nl_mz].append(nl_label)
                    nl_formula = nl_dct[nl_label].get("formula", "")
                    if nl_formula not in raw_residue_nl_dct:
                        raw_residue_nl_dct[nl_formula] = {nl_label: nl_dct[nl_label]}
                    else:
                        raw_residue_nl_dct[nl_formula].update({nl_label: nl_dct[nl_label]})
            else:
                pass
                # print(f'Adduct {pr_adduct} / PR m/z {pr_mz} NOT in range...')
    lipid_space = {}
    lipid_idx_lst = natsorted(list(raw_lipid_space_dct.keys()))
    for lipid_idx in lipid_idx_lst:
        lipid_space[lipid_idx] = raw_lipid_space_dct[lipid_idx]

    residue_class_space = {}
    res_class_idx_lst = natsorted(list(raw_residue_class_dct.keys()))
    for res_class_idx in res_class_idx_lst:
        residue_class_space[res_class_idx] = raw_residue_class_dct[res_class_idx]
    residue_frag_space = {}
    res_frag_idx_lst = natsorted(list(raw_residue_frag_dct.keys()))
    for res_frag_idx in res_frag_idx_lst:
        residue_frag_space[res_frag_idx] = raw_residue_frag_dct[res_frag_idx]
    residue_nl_space = {}
    res_nl_idx_lst = natsorted(list(raw_residue_nl_dct.keys()))
    for res_nl_idx in res_nl_idx_lst:
        residue_nl_space[res_nl_idx] = raw_residue_nl_dct[res_nl_idx]

    del raw_lipid_space_dct
    del raw_residue_class_dct
    del raw_residue_frag_dct
    del raw_residue_nl_dct

    search_space_dct = {
        "lipid_space": lipid_space,
        "residue_space": {
            "class_space": residue_class_space,
            "frag_space": residue_frag_space,
            "nl_space": residue_nl_space,
        },
    }

    return search_space_dct


def build_lipid_search_space(lipid_space: dict, mz_range: list = None):
    if mz_range is None:
        mz_min = 200
        mz_max = 1200
    elif len(mz_range) == 2:
        mz_min = mz_range[0]
        mz_max = mz_range[1]
    else:
        mz_min = 200
        mz_max = 1200

    mz_lipid_space = {}
    reform = {
        (outerKey, innerKey): values
        for outerKey, innerDict in lipid_space.items()
        for innerKey, values in innerDict.items()
    }
    lipid_space_df = pd.DataFrame.from_dict(reform, orient="index")
    mz_lipid_space_df = lipid_space_df.query(f"{mz_min} <= mz <= {mz_max}").copy()
    mz_lipid_idx_lst = mz_lipid_space_df.index.values.tolist()
    elem_lst = list(set([idx[0] for idx in mz_lipid_idx_lst]))
    for elem in elem_lst:
        elem_info = lipid_space.get(elem, {})
        if elem:
            mz_lipid_space[elem] = elem_info

    return mz_lipid_space


def space_to_ranges(
    space: dict,
    ppm: int = 20,
):
    min_ppm_factor = 1 - ppm * 0.000001
    max_ppm_factor = 1 + ppm * 0.000001

    neg_mz_lst = []
    neg_formula_lst = []
    neg_label_lst = []
    pos_mz_lst = []
    pos_formula_lst = []
    pos_label_lst = []
    for formula in space:
        formula_info = space[formula]

        for label in formula_info:
            if label.endswith("-"):
                neg_mz_lst.append(float(formula_info[label].get("mz", 0)))
                neg_formula_lst.append(formula)
                neg_label_lst.append(label)
            elif label.endswith("+"):
                pos_mz_lst.append(float(formula_info[label].get("mz", 0)))
                pos_formula_lst.append(formula)
                pos_label_lst.append(label)
            else:
                pass

    neg_space_df = pd.DataFrame(
        {
            "mz": neg_mz_lst,
            "formula": neg_formula_lst,
            "label": neg_label_lst,
        }
    )
    pos_space_df = pd.DataFrame(
        {
            "mz": pos_mz_lst,
            "formula": pos_formula_lst,
            "label": pos_label_lst,
        }
    )

    neg_space_df["min_mz"] = neg_space_df["mz"] * min_ppm_factor
    neg_space_df["max_mz"] = neg_space_df["mz"] * max_ppm_factor
    pos_space_df["min_mz"] = pos_space_df["mz"] * min_ppm_factor
    pos_space_df["max_mz"] = pos_space_df["mz"] * max_ppm_factor

    neg_range_dct = {
        "min_mz": neg_space_df["min_mz"].min(),
        "max_mz": neg_space_df["max_mz"].max(),
    }
    pos_range_dct = {
        "min_mz": pos_space_df["min_mz"].min(),
        "max_mz": pos_space_df["max_mz"].max(),
    }

    ranges_info = {
        "neg": {"space": neg_space_df, "mz_range": neg_range_dct},
        "pos": {"space": pos_space_df, "mz_range": pos_range_dct},
    }

    return ranges_info


def build_search_ranges(
        space: dict,
        ms1_ppm: int = 20,
        ms2_ppm: int = 20,
):

    search_ranges = {"neg": {}, "pos": {}}

    lipid_space = space.get("lipid_space", None)
    residue_space = space.get("residue_space", None)
    if isinstance(residue_space, dict):
        class_space = residue_space.get("class_space", None)
        frag_space = residue_space.get("frag_space", None)
        nl_space = residue_space.get("nl_space", None)
    else:
        class_space = None
        frag_space = None
        nl_space = None

    if isinstance(lipid_space, dict) and lipid_space:
        lipid_ranges = space_to_ranges(lipid_space, ppm=ms1_ppm)
        search_ranges["neg"]["pr"] = lipid_ranges.get("neg", {})
        search_ranges["pos"]["pr"] = lipid_ranges.get("pos", {})
    if isinstance(class_space, dict) and class_space:
        class_ranges = space_to_ranges(class_space, ppm=ms2_ppm)
        search_ranges["neg"]["class"] = class_ranges.get("neg", {})
        search_ranges["pos"]["class"] = class_ranges.get("pos", {})
    if isinstance(frag_space, dict) and frag_space:
        frag_ranges = space_to_ranges(frag_space, ppm=ms2_ppm)
        search_ranges["neg"]["frag"] = frag_ranges.get("neg", {})
        search_ranges["pos"]["frag"] = frag_ranges.get("pos", {})
    if isinstance(nl_space, dict) and nl_space:
        nl_ranges = space_to_ranges(nl_space, ppm=ms2_ppm)
        search_ranges["neg"]["nl"] = nl_ranges.get("neg", {})
        search_ranges["pos"]["nl"] = nl_ranges.get("pos", {})

    return search_ranges


# def build_residue_ranges(
#     residue_space: dict,
#     ppm: int = 100,
# ):
#     min_ppm_factor = 1 - ppm * 0.000001
#     max_ppm_factor = 1 + ppm * 0.000001
#
#     res_space_dct = {}
#     res_range_dct = {}
#     for sub_res in ["frag_space", "nl_space"]:
#         sub_res_space = residue_space.get(sub_res, {})
#         sub_res_mz_lst = []
#         sub_res_label_lst = []
#         for sub_res_mz_str in sub_res_space:
#             sub_res_mz_lst.append(float(sub_res_mz_str))
#             sub_res_label_lst.append(sub_res_space[sub_res_mz_str])
#
#         sub_res_space_df = pd.DataFrame(
#             {"mz": sub_res_mz_lst, "label": sub_res_label_lst}
#         )
#
#         sub_res_space_df["min_mz"] = sub_res_space_df["mz"] * min_ppm_factor
#         sub_res_space_df["max_mz"] = sub_res_space_df["mz"] * max_ppm_factor
#
#         res_range_dct[sub_res] = {
#             "min_mz": sub_res_space_df["min_mz"].min(),
#             "max_mz": sub_res_space_df["max_mz"].max(),
#         }
#
#         res_space_dct[sub_res] = sub_res_space_df
#
#     res_space_info = {"residue_space": res_space_dct, "residue_range": res_range_dct}
#
#     return res_space_info


if __name__ == "__main__":

    import json
    import time

    # prediction_file = r"../config/usr_epilipidome_PE_PC_lite.json"
    # prediction_file = r"../config/usr_epilipidome_PE_PC.json"
    # prediction_file = r"../config/lite2_epilipidome_PC_max2O_1keto.json"
    # prediction_file = r"../config/lite2_epilipidome_PC_PE_max2O_1keto.json"
    prediction_file = r"../../config/epilipidome_PLs_default_max4O_1keto.json"
    # prediction_file = r"../config/lite3_epilipidome_LPL_PL_TG_CE_max2O_1keto.json"
    usr_mz_range = [600.0, 1000.0]
    t0 = time.time()
    with open(prediction_file, mode="r", encoding="utf-8") as fp:
        usr_epilipidome = json.load(fp)
    print(f"Predicted structures: {len(list(usr_epilipidome.keys()))}")
    t1 = time.time()
    usr_search_space = build_search_space(
        usr_epilipidome, charge_modes=["neg", "pos"],
        adducts=["[M-H]-", "[M+HCOO]-", "[M+Na]+"], mz_range=usr_mz_range
    )
    print(f"Generated search space: {list(usr_search_space.keys())}")
    t2 = time.time()
    with open(
        prediction_file[:-5] + "_search_space.json", mode="w", encoding="utf-8"
    ) as fp:
        json.dump(usr_search_space, fp)
    t3 = time.time()

    print(t1 - t0)
    print(t2 - t1)
    print(t3 - t2)
