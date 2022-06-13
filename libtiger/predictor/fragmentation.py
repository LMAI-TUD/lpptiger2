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

from libtiger.utils.calculator import get_nl_formula, get_exact_mass, get_formula
from libtiger.utils.default_values import (
    charge_cfg_dct,
    specific_ions_dct,
    # ion_mode_elem,
)


def check_pr_adduct(lipid_class: str, pr_name: str, pr_adduct: str):
    """
    Check if the adduct fit to precursor e.g. PC with COOH has [M-H]- instead of [M+HCOO]-
    """
    is_valid_pr_adduct = False
    if lipid_class.endswith("PC"):
        if re.findall(r"COOH", pr_name):
            if pr_adduct == "[M-H]-":
                is_valid_pr_adduct = True
            else:
                pass
        else:
            if pr_adduct in ["[M+HCOO]-", "[M+CH3COO]-"]:
                is_valid_pr_adduct = True
            else:
                pass
    else:
        is_valid_pr_adduct = True

    return is_valid_pr_adduct


def check_res_adduct(lipid_info_dct: dict, lipid_class: str, pr_adduct: str):
    res_lst = []
    for fa in ["FA1", "FA2", "FA3", "FA4"]:
        if fa in lipid_info_dct:
            res_lst.append(fa)
    lipid_class_chg_adduct = charge_cfg_dct.get(lipid_class, {}).get(pr_adduct, [])
    for res in res_lst:
        res_mz_info = lipid_info_dct[res].get("mz_info", {}).copy()
        filtered_res_mz_info = {}
        for res_adduct in res_mz_info:
            if res_adduct in lipid_class_chg_adduct:
                is_valid_adduct = False
                if "-H2O" in res_adduct:
                    oh_rgx = re.compile(r"(.*[<].*)(\d)?([^CO]OH)(.*)")
                    oh_chk = oh_rgx.match(res)
                    if oh_chk and oh_chk.groups()[2] == "OH":
                        is_valid_adduct = True
                else:
                    is_valid_adduct = True
                if is_valid_adduct:
                    filtered_res_mz_info[res_adduct] = res_mz_info[res_adduct]
                else:
                    pass
            else:
                pass
        lipid_info_dct[res]["mz_info"] = filtered_res_mz_info

    return lipid_info_dct


def update_nl_hg(nl_info: dict, pr_formula: str, pr_adduct: str) -> dict:
    neat_nl_formula = nl_info.get("formula")
    pr_nl_formula = get_nl_formula(pr_formula, neat_nl_formula, pr_adduct=pr_adduct, nl_ion_adduct=pr_adduct)
    pr_nl_mz = get_exact_mass(pr_nl_formula)
    nl_info.update({"formula": pr_nl_formula, "mz": pr_nl_mz, "adduct": pr_adduct})

    return nl_info


def update_frag_res(
    lipid_info_dct: dict, lipid_class: str, res_lst: list, ion_adduct_lst: list
):
    frag_dct = {}
    for res in res_lst:
        res_name = lipid_info_dct[res].get("name", "")
        res_ion_mz_info = lipid_info_dct[res].get("mz_info", {}).copy()
        if not res_ion_mz_info:
            print("empty res_ion_mz_info")
        res_ion_mz_info_keys = tuple(res_ion_mz_info.keys())
        adducts_count = 0
        for res_ion_adduct in res_ion_mz_info_keys:
            if res_ion_adduct.startswith('FA'):
                print("!error! res_ion_adduct")
            has_adduct = False
            if res_ion_adduct in ion_adduct_lst:
                if re.match(r'.*H2O.*', res_ion_adduct) and re.match(r'.*OH.*', res_name):
                    has_adduct = True
                elif re.match(r'.*H2O.*', res_ion_adduct) and not re.match(r'.*OH.*', res_name):
                    has_adduct = False
                elif not re.match(r'.*H2O.*', res_ion_adduct) and not re.match(r'.*OH.*', res_name):
                    has_adduct = True
                elif not re.match(r'.*H2O.*', res_ion_adduct) and re.match(r'.*OH.*', res_name):
                    has_adduct = True
                else:
                    pass
                if has_adduct:
                    res_adduct_label = f"{res_name}#{res_ion_adduct}"
                    res_adduct_info = {
                        "name": res_name,
                        "adduct": res_ion_adduct,
                        "label": res_adduct_label,
                        "formula": res_ion_mz_info[res_ion_adduct].get("formula"),
                        "mz": res_ion_mz_info[res_ion_adduct].get("mz"),
                    }
                    # todo zhixu.ni@tu-dresden.de: add smiles
                    res_adduct_info.update(res_ion_mz_info[res_ion_adduct])
                    frag_dct[res_adduct_label] = res_adduct_info
                    adducts_count += 1
                else:
                    # print("skip  ", res_name, res_ion_adduct)
                    pass
            else:
                pass
        if adducts_count < 1:
            print("!error! adducts_count", res_ion_mz_info_keys, res_name, lipid_class)

    return frag_dct


def update_nl_res_pl(
    lipid_info_dct: dict, lipid_class: str, res_lst: list, pr_formula: str, pr_adduct: str
):
    lpl_dct = {}
    # todo zhixu.ni@Uni-leipzig.de: add other lipid classes
    if lipid_class.startswith("P") and len(res_lst) == 2:
        for res in res_lst:
            if res.endswith("1"):
                res_name = lipid_info_dct[f"{res[:-1]}2"].get("name").strip("FA")
            elif res.endswith("2"):
                res_name = lipid_info_dct[f"{res[:-1]}1"].get("name").strip("FA")
            else:
                res_name = None
            if res_name:
                res_info = lipid_info_dct[res]
                res_formula = res_info.get("formula")
                res_no_w_formula = get_nl_formula(res_formula, "H2O", pr_adduct="", nl_ion_adduct="")

                lpl_no_w_formula = get_nl_formula(
                    pr_formula, res_formula, pr_adduct=pr_adduct, nl_ion_adduct="[M-H]-"
                )
                lpl_formula = get_nl_formula(
                    pr_formula, res_no_w_formula, pr_adduct=pr_adduct, nl_ion_adduct="[M-H]-"
                )
                lpl_mz = get_exact_mass(lpl_formula)
                lpl_no_w_mz = get_exact_mass(lpl_no_w_formula)
                if lipid_class.startswith("L"):
                    pass
                else:
                    lipid_class = "L" + lipid_class

                lpl_name = f"{lipid_class}({res_name})"
                res_lpl_dct = {
                    f"{lpl_name}#[M-H]-": {
                        "name": lpl_name,
                        "adduct": "[M-H]-",
                        "label": f"{lpl_name}#[M-H]-",
                        "formula": lpl_formula,
                        "mz": lpl_mz,
                    },
                    f"{lpl_name}#[M-H2O-H]-": {
                        "name": lpl_name,
                        "adduct": "[M-H2O-H]-",
                        "label": f"{lpl_name}#[M-H2O-H]-",
                        "formula": lpl_no_w_formula,
                        "mz": lpl_no_w_mz,
                    },
                }
                if lipid_class in ["PC", "LPC"]:
                    lpc_no_w_formula = get_nl_formula(
                        pr_formula, res_formula, pr_adduct=pr_adduct, nl_ion_adduct="[M-CH3]-"
                    )
                    lpc_formula = get_nl_formula(
                        pr_formula, res_no_w_formula, pr_adduct=pr_adduct, nl_ion_adduct="[M-CH3]-"
                    )
                    lpc_mz = get_exact_mass(lpc_formula)
                    lpc_no_w_mz = get_exact_mass(lpc_no_w_formula)
                    res_lpl_dct.update(
                        {
                            f"{lpl_name}#[M-CH3]-": {
                                "name": lpl_name,
                                "adduct": "[M-CH3]-",
                                "label": f"{lpl_name}#[M-CH3]-",
                                "formula": lpc_formula,
                                "mz": lpc_mz,
                            },
                            f"{lpl_name}#[M-CH3-H2O]-": {
                                "name": lpl_name,
                                "adduct": "[M-CH3-H2O]-",
                                "label": f"{lpl_name}#[M-CH3-H2O]-",
                                "formula": lpc_no_w_formula,
                                "mz": lpc_no_w_mz,
                            },
                        }
                    )
                lpl_dct.update(res_lpl_dct)

    return lpl_dct


def update_nl_res_tg(
    lipid_info_dct: dict, nl_lipid_class: str, res_lst: list, pr_formula: str, pr_adduct: str
):
    dg_dct = {}
    res_names_lst = [lipid_info_dct[r].get("name").strip("FA") for r in res_lst]
    if nl_lipid_class.startswith("TG") and len(res_lst) == 3:
        for res in res_lst:
            res_name = lipid_info_dct[res].get("name").strip("FA")
            other_res_name_lst = [ot_name for ot_name in res_names_lst]  # safe copy of list
            other_res_name_lst.remove(res_name)
            nl_res_name_lst = natsorted(other_res_name_lst)

            if res_name:
                res_info = lipid_info_dct[res]
                res_formula = res_info.get("formula")
                res_no_w_formula = get_nl_formula(res_formula, "H2O", pr_adduct="", nl_ion_adduct="")

                dg_no_w_formula = get_nl_formula(
                    pr_formula, res_formula,  pr_adduct=pr_adduct, nl_ion_adduct="[M+H]+"
                )
                dg_formula = get_nl_formula(
                    pr_formula, res_no_w_formula, pr_adduct=pr_adduct, nl_ion_adduct="[M+H]+"
                )
                dg_no_w_sod_formula = get_nl_formula(
                    pr_formula, res_formula, pr_adduct=pr_adduct, nl_ion_adduct="[M+Na]+"
                )
                dg_sod_formula = get_nl_formula(
                    pr_formula, res_no_w_formula, pr_adduct=pr_adduct, nl_ion_adduct="[M+Na]+"
                )
                dg_mz = get_exact_mass(dg_formula)
                dg_no_w_mz = get_exact_mass(dg_no_w_formula)
                dg_sod_mz = get_exact_mass(dg_sod_formula)
                dg_no_w_sod_mz = get_exact_mass(dg_no_w_sod_formula)
                nl_lipid_class = "DG"

                dg_name = f"{nl_lipid_class}({'_'.join(nl_res_name_lst)})"
                res_lpl_dct = {
                    f"{dg_name}#[M+H]+": {
                        "name": dg_name,
                        "adduct": "[M+H]+",
                        "label": f"{dg_name}#[M+H]+",
                        "formula": dg_formula,
                        "mz": dg_mz,
                    },
                    f"{dg_name}#[M-H2O+H]+": {
                        "name": dg_name,
                        "adduct": "[M-H2O+H]+",
                        "label": f"{dg_name}#[M-H2O+H]+",
                        "formula": dg_no_w_formula,
                        "mz": dg_no_w_mz,
                    },
                    f"{dg_name}#[M+Na]+": {
                        "name": dg_name,
                        "adduct": "[M+Na]+",
                        "label": f"{dg_name}#[M+Na]+",
                        "formula": dg_sod_formula,
                        "mz": dg_sod_mz,
                    },
                    f"{dg_name}#[M-H2O+Na]+": {
                        "name": dg_name,
                        "adduct": "[M-H2O+Na]+",
                        "label": f"{dg_name}#[M-H2O+Na]+",
                        "formula": dg_no_w_sod_formula,
                        "mz": dg_no_w_sod_mz,
                    },
                }

                dg_dct.update(res_lpl_dct)

    return dg_dct


def get_fingerprint(lipid_dct: dict) -> list:
    pr_name = lipid_dct.get("name", "")
    pr_mz = lipid_dct.get("mz", 0.0)
    ions = lipid_dct.get("ions", {})
    fp_lst = [pr_mz]
    fp_oh_lst = []
    fp_cooh_lst = []
    oh_rgx = re.compile(r"(.*[<][^CO]*)(\dOH|OH)(.*)")
    ooh_rgx = re.compile(r"(.*[<][^C]*)(\dOOH|OOH)(.*)")
    for ion_typ in ions:
        sub_gp_ions = ions[ion_typ]
        if ion_typ == "class_info":
            for ion_label in sub_gp_ions:
                ion_mz = sub_gp_ions[ion_label].get("mz")
                fp_lst.append(ion_mz)
                oh_chk = oh_rgx.match(pr_name)
                if oh_chk and oh_chk.groups()[1] is not None:
                    oh_seg = oh_chk.groups()[1]
                    if oh_seg[0] == "O":
                        oh_count = 1
                    elif re.match(r"\d", oh_seg[0]):
                        try:
                            oh_count = int(oh_seg[0])
                        except TypeError:
                            oh_count = 0
                    else:
                        oh_count = 0
                    if oh_count > 0:
                        for oh_c in range(1, oh_count + 1):
                            fp_oh_lst.append(ion_mz - oh_c * 18.010565)  # -H2O
                    else:
                        pass
            else:
                for ion_label in sub_gp_ions:
                    ion_mz = sub_gp_ions[ion_label].get("mz")
                    fp_lst.append(ion_mz)
        else:
            for ion_label in sub_gp_ions:
                ion_mz = sub_gp_ions[ion_label].get("mz")
                ion_name = sub_gp_ions[ion_label].get("name")
                fp_oh_lst.append(ion_mz)
                oh_chk = oh_rgx.match(ion_name)
                if oh_chk and oh_chk.groups()[1] is not None:
                    oh_seg = oh_chk.groups()[1]
                    if oh_seg[0] == "O":
                        oh_count = 1
                    elif re.match(r"\d", oh_seg[0]):
                        try:
                            oh_count = int(oh_seg[0])
                        except TypeError:
                            oh_count = 0
                    else:
                        oh_count = 0
                    if oh_count > 0:
                        for oh_c in range(1, oh_count + 1):
                            fp_oh_lst.append(ion_mz - oh_c * 18.010565)  # -H2O
                    else:
                        pass
                else:
                    pass
                ooh_chk = ooh_rgx.match(ion_name)
                if ooh_chk and ooh_chk.groups()[1] is not None:
                    ooh_seg = ooh_chk.groups()[1]
                    if ooh_seg[0] == "O":
                        ooh_count = 1
                    elif re.match(r"\d", ooh_seg[0]):
                        try:
                            ooh_count = int(ooh_seg[0])
                        except TypeError:
                            ooh_count = 0
                    else:
                        ooh_count = 0
                    if ooh_count > 0:
                        for ooh_c in range(1, ooh_count + 1):
                            fp_oh_lst.append(ion_mz - ooh_c * 18.010565)  # -H2O
                    else:
                        pass
                else:
                    pass

                if "COOH" in ion_name:
                    for fp_oh_mz in fp_oh_lst:
                        fp_cooh_lst.append(fp_oh_mz - 43.989830)  # -CO2
                else:
                    pass

                fp_lst.extend(fp_oh_lst)
                fp_lst.extend(fp_cooh_lst)
    # add other peaks in msp
    msp_mz_lst = lipid_dct.get("msp", {}).get("mz", [])
    fp_lst.extend(msp_mz_lst)
    fp_lst = sorted(list(set(fp_lst)))
    fp_lst = [round(fp_val, 6) for fp_val in fp_lst if fp_val > 0]

    return fp_lst


def get_matched_peak(
    peak_label: str,
    raw_peak_label: str,
    sum_ions: dict,
    msp_frag_pattern_info: dict,
    msp_dct: dict,
):
    if peak_label in sum_ions:
        msp_dct["mz"].append(sum_ions[peak_label].get("mz", 0.0))
        msp_dct["i"].append(msp_frag_pattern_info.get(raw_peak_label))
        msp_dct["label"].append(peak_label)
        msp_dct["remark"].append(
            f"{peak_label}#{sum_ions[peak_label].get('formula', '')}"
        )
    else:
        peak_base = peak_label.split("#")[0]
        peak_chg = peak_label.split("#")[1]
        base_label = f"{peak_base}#[M-H]-"
        base_formula = sum_ions.get(base_label, {}).get("formula", "")
        mz = get_exact_mass(base_formula.strip("+-") + "H", charge=peak_chg)
        formula = get_formula(base_formula.strip("+-") + "H", charge=peak_chg)
        if mz and formula:
            msp_dct["mz"].append(mz)
            msp_dct["i"].append(msp_frag_pattern_info.get(raw_peak_label))
            msp_dct["label"].append(peak_label)
            msp_dct["remark"].append(f"{peak_label}#{formula}")
        else:
            # print("Cannot find matched MSP peak config.")
            pass

    return msp_dct


def get_msp(
    lipid_dct: dict,
    class_info: dict,
    frag_patterns: dict,
    lipid_class: str,
    pr_adduct: str,
):
    msp_dct = {}
    lipid_name = lipid_dct.get("name", "")
    mod_typ = lipid_dct.get("mod_type", "")
    pr_mz = lipid_dct.get("mz", 0.0)
    pr_formula = lipid_dct.get("formula", "")
    lipid_info = lipid_dct.get("info", "")
    fa1_name = lipid_info.get("FA1", {}).get("name", "")
    fa2_name = lipid_info.get("FA2", {}).get("name", "")
    fa3_name = lipid_info.get("FA3", {}).get("name", "")
    msp_frag_pattern_info = {}
    has_matched_pattern = False
    if lipid_class in frag_patterns:
        class_frag_patterns_dct = frag_patterns[lipid_class]
        pr_adduct_patterns_dct = class_frag_patterns_dct.get(pr_adduct, {})
        pr_mod_patterns_dct = pr_adduct_patterns_dct.get(mod_typ, {})

        for pr_rgx_str in pr_mod_patterns_dct:
            if re.match(pr_rgx_str, lipid_name):
                msp_frag_pattern_info = pr_mod_patterns_dct.get(pr_rgx_str, {})
                has_matched_pattern = True

        if has_matched_pattern:
            pass
        else:
            pr_mod_patterns_dct = pr_adduct_patterns_dct.get("unmod", {})
            umod_rgx_str_lst = list(pr_mod_patterns_dct.keys())
            if umod_rgx_str_lst:
                msp_frag_pattern_info = pr_mod_patterns_dct.get(umod_rgx_str_lst[0], "")
            else:
                pass

    if msp_frag_pattern_info:
        msp_dct = {
            "mz": [],
            "i": [],
            "label": [],
            "remark": [],
        }
        ions_dct = lipid_dct.get("ions", {})
        sum_ions_dct = {}
        for sub_group in ions_dct:
            sum_ions_dct.update(ions_dct.get(sub_group))
        for raw_peak_label in msp_frag_pattern_info:
            if raw_peak_label == f"PR#{pr_adduct}":
                msp_dct["mz"].append(pr_mz)
                msp_dct["i"].append(msp_frag_pattern_info.get(raw_peak_label))
                msp_dct["label"].append(raw_peak_label)
                msp_dct["remark"].append(f"{lipid_name}#{pr_adduct}|{pr_formula}")
            elif re.match(r"^FA1.*$", raw_peak_label):
                peak_label = re.sub(r"FA1", fa1_name, raw_peak_label)
                msp_dct = get_matched_peak(
                    peak_label,
                    raw_peak_label,
                    sum_ions_dct,
                    msp_frag_pattern_info,
                    msp_dct,
                )
            elif re.match(r".*\(FA1.*$", raw_peak_label):
                peak_label = re.sub(r"FA1", fa2_name.strip("FA"), raw_peak_label)
                msp_dct = get_matched_peak(
                    peak_label,
                    raw_peak_label,
                    sum_ions_dct,
                    msp_frag_pattern_info,
                    msp_dct,
                )
            elif re.match(r"^FA2.*$", raw_peak_label):
                peak_label = re.sub(r"FA2", fa2_name, raw_peak_label)
                msp_dct = get_matched_peak(
                    peak_label,
                    raw_peak_label,
                    sum_ions_dct,
                    msp_frag_pattern_info,
                    msp_dct,
                )
            elif re.match(r".*\(FA2.*$", raw_peak_label):
                peak_label = re.sub(r"FA2", fa1_name.strip("FA"), raw_peak_label)
                msp_dct = get_matched_peak(
                    peak_label,
                    raw_peak_label,
                    sum_ions_dct,
                    msp_frag_pattern_info,
                    msp_dct,
                )
            elif re.match(r"^FA3.*$", raw_peak_label):
                peak_label = re.sub(r"FA3", fa2_name, raw_peak_label)
                msp_dct = get_matched_peak(
                    peak_label,
                    raw_peak_label,
                    sum_ions_dct,
                    msp_frag_pattern_info,
                    msp_dct,
                )
            elif re.match(r".*\(FA3.*$", raw_peak_label):
                peak_label = re.sub(r"FA3", fa1_name.strip("FA"), raw_peak_label)
                msp_dct = get_matched_peak(
                    peak_label,
                    raw_peak_label,
                    sum_ions_dct,
                    msp_frag_pattern_info,
                    msp_dct,
                )
            elif raw_peak_label in class_info:
                msp_dct = get_matched_peak(
                    raw_peak_label,
                    raw_peak_label,
                    sum_ions_dct,
                    msp_frag_pattern_info,
                    msp_dct,
                )
            else:
                pass

    return msp_dct


def generate_fragments(lipid_dct: dict, frag_patterns: dict):

    pr_dct = {}

    pr_info_dct = lipid_dct.copy()

    pr_name = lipid_dct.get("name", "")
    pr_neutral_formula = lipid_dct.get("formula", "")
    lipid_info_dct = lipid_dct.get("info", {}).copy()
    lipid_class = lipid_info_dct.get("lipid_class", "")
    pr_adducts_dct = charge_cfg_dct.get(lipid_class, {})
    pr_adducts = list(pr_adducts_dct.keys())

    hg_info = specific_ions_dct.get(lipid_class, {})
    hg_frag_info = hg_info.get("frag_info", {})
    hg_nl_info = hg_info.get("nl_info", {})

    for pr_adduct in pr_adducts:

        is_valid_pr_adduct = check_pr_adduct(lipid_class, pr_name, pr_adduct)

        if is_valid_pr_adduct:
            class_info = hg_frag_info
            for hg_nl in hg_nl_info:
                class_info[hg_nl] = update_nl_hg(
                    hg_nl_info[hg_nl].copy(), pr_neutral_formula, pr_adduct
                )

            pr_mz_info = lipid_dct.get("mz_info", {}).get(pr_adduct)
            pr_ion_formula = pr_mz_info.get("formula")
            ion_adduct_lst = pr_adducts_dct.get(pr_adduct, [])
            res_lst = []
            for fa in ["FA1", "FA2", "FA3"]:
                if fa in lipid_info_dct:
                    res_lst.append(fa)

            frag_info = update_frag_res(lipid_info_dct, lipid_class, res_lst, ion_adduct_lst)
            if lipid_class.startswith("P"):
                nl_info = update_nl_res_pl(lipid_info_dct, lipid_class, res_lst, pr_ion_formula, pr_adduct)
            elif lipid_class.startswith("TG"):
                nl_info = update_nl_res_tg(lipid_info_dct, lipid_class, res_lst, pr_ion_formula, pr_adduct)
            else:
                nl_info = {}
            pr_info_dct["adduct"] = pr_adduct
            pr_info_dct["formula"] = pr_ion_formula
            pr_info_dct["mz"] = pr_mz_info.get("mz")
            pr_info_dct["info"] = check_res_adduct(
                lipid_info_dct, lipid_class, pr_adduct
            )

            # if frag_info:
            #     # print("   ", pr_name, pr_adduct)
            #     pass
            # else:
            #     print("BAD", pr_name, pr_adduct)
            #     # print(pr_name)

            rank_ions = {
                "class_info": class_info,
                "frag_info": frag_info,
                "nl_info": nl_info,
            }

            if lipid_class == "CE":
                ce_frag = rank_ions.get("frag_info", {}).copy()
                ce_frag.update(class_info)
                rank_ions["frag_info"] = ce_frag

            pr_info_dct["ions"] = rank_ions
            pr_info_dct["msp"] = get_msp(
                pr_info_dct,
                class_info,
                frag_patterns,
                lipid_class=lipid_class,
                pr_adduct=pr_adduct,
            )
            pr_info_dct["fingerprints"] = get_fingerprint(pr_info_dct)
            if "mz_info" in pr_info_dct:
                del pr_info_dct["mz_info"]
            if "exact_mass" in pr_info_dct:
                del pr_info_dct["exact_mass"]
            pr_info_dct["neutral_formula"] = (
                lipid_dct.get("mz_info", {}).get("M", {}).get("formula")
            )
            pr_info_dct["neutral_mass"] = (
                lipid_dct.get("mz_info", {}).get("M", {}).get("mz")
            )
            pr_dct = {f"{pr_name}#{pr_adduct}": pr_info_dct}

    return pr_dct


if __name__ == "__main__":

    import json

    from libtiger.reader.config_parser import load_frag_pattern

    # prediction_file = r"../config/temp.json"
    # prediction_file = r"../config/usr_epilipidome_PLs_lite.json"
    # prediction_file = r"../config/usr_epilipidome_PE_PC_lite.json"
    prediction_file = r"../../config/usr_epilipidome_PE_lite.json"
    usr_frag_patterns_xlsx = r"../config/4_FragmentationPatterns.xlsx"
    usr_frag_patterns = load_frag_pattern(usr_frag_patterns_xlsx)
    with open(prediction_file, mode="r", encoding="utf-8") as fp:
        usr_epilipidome = json.load(fp)
    print(len(list(usr_epilipidome.keys())))
    for usr_lipid in usr_epilipidome:
        usr_fragments_info = generate_fragments(
            usr_epilipidome[usr_lipid], frag_patterns=usr_frag_patterns
        )
        print(usr_fragments_info)
    print("FIN")
