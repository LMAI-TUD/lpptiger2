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
from typing import List

import pandas as pd

from libtiger.utils.calculator import calc_smi


fa_rgx = re.compile(r"(?P<C>\d{1,2})(:)(?P<DB>\d)")
temp_file_path = r"../temp/predicted_"


def load_modifications(file: str, mod_level: int = 1):

    mod_df = pd.read_csv(file, index_col=0).T
    mod_df = mod_df.astype(
        {
            "OXIDATIONLEVEL": "int32",
            "OAP": "int32",
            "OCP": "int32",
            "DB": "int32",
            "OH": "int32",
            "KETO": "int32",
            "CHO": "int32",
            "COOH": "int32",
            "EPOXY": "int32",
            "OOH": "int32",
        }
    )
    mod_selected_df = mod_df[mod_df["OXIDATIONLEVEL"] == mod_level]
    # mod_dct = mod_selected_df.to_dict(orient="index")
    # print(mod_dct)

    return mod_df


def load_residues(file: str):
    sum_res_dct = {}
    res_df = pd.read_csv(file)
    res_df.drop_duplicates("Residues", inplace=True)
    res_df.set_index("Residues", drop=True, inplace=True)
    smi_key = None
    df_col_lst = res_df.columns.values.tolist()
    for k in ["smi", "SMILES", "smiles"]:
        if k in df_col_lst:
            smi_key = k
    res_df.rename(columns={smi_key: "smiles"}, inplace=True)
    if smi_key:
        res_dct = res_df.to_dict(orient="index")
        for res in res_dct:
            res_smi = res_dct[res].get("smiles", "")
            tmp_res_dct = {"name": res, "mod_type": "unmod", "smiles": res_smi}
            tmp_res_dct.update(calc_smi(lipid_class="FA", smiles=res_smi))
            tmp_res_dct["origins"] = [res]
            sum_res_dct[res] = tmp_res_dct
    else:
        sum_res_dct = {}

    return sum_res_dct


def load_frag_pattern(file: str):

    xls = pd.ExcelFile(file)
    sheet_names = xls.sheet_names
    frag_pattern_info = {}
    for lipid_class in ["PA", "PC", "PE", "PG", "PS"]:
    # for lipid_class in ["PA", "PC", "PE", "PG", "PS", "TG", "CE"]:
        if lipid_class in sheet_names:
            frag_pattern_df = pd.read_excel(file, sheet_name=lipid_class)
            frag_pattern_df.set_index("FA_pattern", drop=True, inplace=True)
            chg_type_patterns = {}
            chg_lst = list(set(frag_pattern_df["PR_CHARGE"].values.tolist()))
            for chg_typ in chg_lst:
                chg_df = frag_pattern_df[frag_pattern_df["PR_CHARGE"] == chg_typ].copy()
                mod_type_patterns = {}
                for mod_typ in ["oap", "ocp", "unmod"]:
                    mod_df = chg_df[chg_df["TYPE"] == mod_typ].copy()
                    raw_frag_pattern_dct = mod_df.to_dict(orient="index")
                    frag_pattern_dct = {}
                    for i_k in raw_frag_pattern_dct:
                        tmp_pattern = {}
                        raw_frag_pattern = raw_frag_pattern_dct[i_k]
                        for frag_k in raw_frag_pattern:
                            if "#" in frag_k:
                                try:
                                    frag_i = int(raw_frag_pattern[frag_k])
                                except TypeError:
                                    frag_i = 0
                                if frag_i > 0:
                                    tmp_pattern[frag_k] = int(raw_frag_pattern[frag_k])
                        rgx_str = f"{lipid_class}\(\d\d?:\d_{i_k}\)"
                        frag_pattern_dct[rgx_str] = tmp_pattern
                    mod_type_patterns[mod_typ] = frag_pattern_dct
                chg_type_patterns[chg_typ] = mod_type_patterns
            frag_pattern_info[lipid_class] = chg_type_patterns

    return frag_pattern_info


# def construct_params(observed_classes: list, cfg: dict, precursors: list):
#     sum_params = {}
#     settings = {
#         "fa_lst_path": r"../config/2_FA_list.csv",
#         "mod_lst_path": r"../config/2_DB_Mod_cfg.csv",
#         "ox_level": 1,
#         "oap_mode": 1,
#         "ocp_mode": 1,
#         "lyso_oap_mode": 1,
#         "lyso_ocp_mode": 1,
#         "prostane_mode": 0,
#         "ox_prostane_mode": 0,
#         "msp_mode": 0,
#         "frag_pattern_path": r"../ConfigurationFiles/5_Fragmentation_Patterns.xlsx",
#         "pl_hg_path": r"../ConfigurationFiles/6_PL_Specific_Signals_cfg.xlsx",
#         "prostane_mod_path": r"../ConfigurationFiles/3_Prostanes_Mod_cfg.csv",
#         "prostane_abbr_path": r"../ConfigurationFiles/4_Prostanes_Abbr_cfg.csv",
#     }
#     shared_params = {
#         "lipid_lst_path": r"{Lipid_list to be oxidized}",
#         "max_oh": cfg.get("max_OH", 3),
#         "max_oxo": cfg.get("max_oxo", 1),
#         "max_ooh": cfg.get("max_OOH", 1),
#         "max_epoxy": cfg.get("max_epoxy", 0),
#     }
#
#     for lipid_class in observed_classes:
#         lipid_class_params = {}
#         lipid_class_params.update(settings)
#         lipid_class_params.update(shared_params)
#         lipid_class_params["lipid_class"] = lipid_class
#
#         precursor_df = precursors[precursors["lipid_class"] == lipid_class]
#         precursor_lst = precursor_df["Converted_Name"].values.tolist()
#
#         lipid_class_params["lipid_lst"] = precursor_lst
#         lipid_class_params["sdf_path"] = f"{temp_file_path}{lipid_class}.sdf"
#         lipid_class_params["msp_path"] = f"{temp_file_path}{lipid_class}.msp"
#
#         sum_params[lipid_class] = lipid_class_params
#
#     return sum_params


# def filter_residues(names: List[str]) -> (dict, dict):
#     """
#     Split saturated FA and unsaturated FA
#     """
#
#     saturated_fa_dct = {}
#     unsaturated_fa_dct = {}
#
#     for res_name in names:
#         fa_match = fa_rgx.match(res_name)
#         if fa_match:
#             fa_groups = fa_match.groupdict()
#             print(fa_groups)
#             if int(fa_groups["DB"]) == 0:
#                 saturated_fa_dct[res_name] = {
#                     "C": int(fa_groups["C"]),
#                     "DB": int(fa_groups["DB"]),
#                 }
#             else:
#                 unsaturated_fa_dct[res_name] = {
#                     "C": int(fa_groups["C"]),
#                     "DB": int(fa_groups["DB"]),
#                 }
#
#     return saturated_fa_dct, unsaturated_fa_dct


if __name__ == "__main__":

    usr_res_file = r"../../config/1_FA_list.csv"
    load_residues(usr_res_file)

    usr_mod_file = r"../../config/2_DB_Mod_cfg.csv"
    load_modifications(usr_mod_file)

    usr_frag_pattern_file = r"../../config/4_Fragmentation_Patterns.xlsx"
    load_frag_pattern(usr_frag_pattern_file)
