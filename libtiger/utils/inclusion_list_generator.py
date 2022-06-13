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
import json
import re

import pandas as pd


def gen_incl(space, origins, export_file):
    headers = [
        "lipid_name",
        "neutral_formula",
        "neutral_mass",
        "charge",
        "charged_formula",
        "mz",
        "oxFA_name",
        "oxFA_neutral_formula",
        "oxFA_neutral_mass",
        "oxFA_charge",
        "oxFA_charged_formula",
        "oxFA_mz",
    ]

    # space = r"/home/ni/sysmedos/tiger2star/libtiger/config/lite_epilipidome_LPL_max2O_1keto.json"
    # space = r"D:\SysMedOs\tigerstar\tiger2star\libtiger\config\lite3_epilipidome_LPL_PL_TG_CE_max2O_1keto.json"

    origin_df = pd.read_excel(origins, engine="openpyxl")
    lipid_col_name = "lipid"
    headers_lst = origin_df.columns.tolist()
    for col in headers_lst:
        if re.match(r'.*lipid.*', col, re.IGNORECASE):
            lipid_col_name = col
    origin_lst = origin_df[lipid_col_name].tolist()
    print(origin_lst)
    found_ox_lst = []

    with open(space, "r") as f:
        space_dct = json.load(f)

    print(len(list(space_dct.keys())))
    print("Space loaded...")

    df_dct = {
        "lipid_name": [],
        "neutral_formula": [],
        "neutral_mass": [],
        "charge": [],
        "charged_formula": [],
        "mz": [],
        "oxFA_name": [],
        "oxFA_neutral_formula": [],
        "oxFA_neutral_mass": [],
        # "oxFA_charge": [],
        # "oxFA_charged_formula": [],
        # "oxFA_mz": [],
    }

    for k in space_dct:
        for k_chg in space_dct[k]:
            k_info = space_dct[k].get(k_chg, {})
            mod_type = k_info.get("mod_type", "")
            is_candidate = 0
            if mod_type and mod_type != "unmod":
                k_origin_lst = k_info.get("origins", [])
                for o in origin_lst:
                    if o in k_origin_lst:
                        is_candidate = 1
                        found_ox_lst.append(o)
                    else:
                        pass
            if is_candidate:
                print(k_info.get("name", ""))
                res_info = k_info.get("info", {})
                res_mod_count = 0
                k_mod_res = ""
                for k_res in res_info:
                    if k_res in ["FA1", "FA2", "FA3"]:
                        res_mod_type = res_info[k_res].get("mod_type", "")
                        if res_mod_type and res_mod_type != "unmod":
                            tmp_mod_res_info = res_info.get(k_res, {})
                            ox_name = tmp_mod_res_info.get("name", "")
                            if "<" in ox_name:
                                res_mod_count += 1
                                k_mod_res = k_res
                if res_mod_count == 1:
                    print(k, 1)
                if res_mod_count == 1 and k_mod_res:
                    mod_res_info = res_info.get(k_mod_res, {})
                    print(k_info.get("name", ""))
                    df_dct["lipid_name"].append(k_info.get("name", ""))
                    print(k_info.get("neutral_formula", ""))
                    df_dct["neutral_formula"].append(k_info.get("neutral_formula", ""))
                    print(k_info.get("neutral_mass", 0))
                    df_dct["neutral_mass"].append(k_info.get("neutral_mass", 0))
                    print(k_info.get("adduct", ""))
                    df_dct["charge"].append(k_info.get("adduct", ""))
                    print(k_info.get("formula", ""))
                    df_dct["charged_formula"].append(k_info.get("formula", ""))
                    print(k_info.get("mz", 0))
                    df_dct["mz"].append(k_info.get("mz", 0))
                    print(mod_res_info.get("name", ""))
                    df_dct["oxFA_name"].append(mod_res_info.get("name", ""))
                    print(mod_res_info.get("formula", ""))
                    df_dct["oxFA_neutral_formula"].append(
                        mod_res_info.get("formula", "")
                    )
                    print(mod_res_info.get("exact_mass", ""))
                    df_dct["oxFA_neutral_mass"].append(
                        mod_res_info.get("exact_mass", "")
                    )
                    res_mz_info = mod_res_info.get("mz_info", "")
                    # for chg in ["[M-H]-", "[M+Na]+"]:
                    #     if chg in res_mz_info:
                    #         df_dct["oxFA_charge"].append(chg)
                    #         df_dct["oxFA_charged_formula"].append(res_mz_info[chg].get("formula", ""))
                    #         df_dct["oxFA_mz"].append(res_mz_info[chg].get("mz", 0))
                    #         break
                    print(
                        k_info.get("name", ""),
                        k_info.get("neutral_mass", 0),
                        mod_res_info.get("name", ""),
                    )
    obs_pr_lst = list(set(found_ox_lst))
    skipped_pr_lst = [pr for pr in origin_lst if pr not in obs_pr_lst]
    df = pd.DataFrame(df_dct)
    print(df.head())
    print("skipped_pr_lst")
    print(skipped_pr_lst)
    df.to_csv(export_file)

    print("FIN")
    return skipped_pr_lst


if __name__ == "__main__":
    usr_space = r"../config/lite4_epilipidome_LPL_PL_TG_CE_max2O_1keto.json"

    in_lst = [
        # r"/home/ni/sysmedos/tiger2star/temp/Palina/cortex/cortex_CE.xlsx",
        # r"/home/ni/sysmedos/tiger2star/temp/Palina/cortex/cortex_GL.xlsx",
        # r"/home/ni/sysmedos/tiger2star/temp/Palina/cortex/cortex_LPL.xlsx",
        # r"/home/ni/sysmedos/tiger2star/temp/Palina/cortex/cortex_PL.xlsx",
        # r"/home/ni/sysmedos/tiger2star/temp/Palina/cereb/cereb_CE.xlsx",
        # r"/home/ni/sysmedos/tiger2star/temp/Palina/cereb/cereb_GL.xlsx",
        # r"/home/ni/sysmedos/tiger2star/temp/Palina/cereb/cereb_LPL.xlsx",
        # r"/home/ni/sysmedos/tiger2star/temp/Palina/cereb/cereb_PL.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\plasma\plasma_neg_LPL.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\plasma\plasma_neg_PL.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\plasma\plasma_pos_CE.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\plasma\plasma_pos_TG.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\cortex\cortex_neg_LPL.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\cortex\cortex_neg_PL.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\cortex\cortex_pos_CE.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\cortex\cortex_pos_TG.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\cereb\cereb_neg_LPL.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\cereb\cereb_neg_PL.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\cereb\cereb_pos_CE.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\cereb\cereb_pos_TG.xlsx",
        r"D:\TUD\Nextcloud\MF\predictions121221\predictions121221_input.xlsx",
    ]
    out_lst = [
        # r"/home/ni/sysmedos/tiger2star/temp/Palina/cortex/ox_cortex_pos_CE.xlsx",
        # r"/home/ni/sysmedos/tiger2star/temp/Palina/cortex/ox_cortex_pos_GL.xlsx",
        # r"/home/ni/sysmedos/tiger2star/temp/Palina/cortex/ox_cortex_neg_LPL.xlsx",
        # r"/home/ni/sysmedos/tiger2star/temp/Palina/cortex/ox_cortex_neg_PL.xlsx",
        # r"/home/ni/sysmedos/tiger2star/temp/Palina/cereb/ox_cereb_pos_CE.xlsx",
        # r"/home/ni/sysmedos/tiger2star/temp/Palina/cereb/ox_cereb_pos_GL.xlsx",
        # r"/home/ni/sysmedos/tiger2star/temp/Palina/cereb/ox_cereb_neg_LPL.xlsx",
        # r"/home/ni/sysmedos/tiger2star/temp/Palina/cereb/ox_cereb_neg_PL.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\plasma\plasma_ox_neg_LPL.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\plasma\plasma_ox_neg_PL.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\plasma\plasma_ox_pos_CE.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\plasma\plasma_ox_pos_TG.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\cortex\cortex_ox_neg_LPL.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\cortex\cortex_ox_neg_PL.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\cortex\cortex_ox_pos_CE.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\cortex\cortex_ox_pos_TG.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\cereb\cereb_ox_neg_LPL.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\cereb\cereb_ox_neg_PL.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\cereb\cereb_ox_pos_CE.xlsx",
        # r"D:\SysMedOs\tigerstar\tiger2star\temp\cereb\cereb_ox_pos_TG.xlsx",
        r"D:\TUD\Nextcloud\MF\predictions121221\predictions121221_output.xlsx",
    ]

    pairs = list(zip(in_lst, out_lst))
    skip_dct = {}
    for p in pairs:
        skip_lst = gen_incl(usr_space, p[0], p[1])
        skip_dct[p[0]] = skip_lst

    print(skip_dct)

    print('FIN')
