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
import itertools
import json
import math
import multiprocessing
import os.path
from multiprocessing import Pool
import time

from natsort import natsorted
import pandas as pd
import regex as re

from libtiger.predictor.fragmentation import generate_fragments
from libtiger.utils.calculator import calc_smi
from libtiger.utils.default_values import pl_hg_dct, res_rgx, unit_rgx
from libtiger.reader.config_parser import load_frag_pattern, load_modifications, load_residues
from libtiger.utils.refine_epilipidome import refine_epilipidome


def fit_mod_count(smi: str, max_mod_cfg: dict = None):

    is_allowed = False

    oh_rgx = re.compile(r"C\(O\)")
    oxo_rgx = re.compile(r"C\(=O\)")
    ooh_rgx = re.compile(r"C\(OO\)")
    epoxy_rgx = re.compile(r"C1C\(O1\)C")
    cho_rgx = re.compile(r"C=O")
    cooh_rgx = re.compile(r"C\(O\)=O")

    act_oh_count = len(oh_rgx.findall(smi))
    act_oxo_count = len(oxo_rgx.findall(smi))
    act_ooh_count = len(ooh_rgx.findall(smi))
    act_epoxy_count = len(epoxy_rgx.findall(smi))
    act_cho_count = len(cho_rgx.findall(smi))
    act_cooh_count = len(cooh_rgx.findall(smi))

    if not max_mod_cfg:
        max_mod_cfg = {
            "max_mod": 2,
            "max_o": 3,
            "max_oh": 2,
            "max_oxo": 1,
            "max_ooh": 1,
            "max_ep": 0,
        }
    max_mod_count = max_mod_cfg.get("max_mod", 2)
    max_o_count = max_mod_cfg.get("max_o", 3)
    max_oh_count = max_mod_cfg.get("max_oh", 2)
    max_oxo_count = max_mod_cfg.get("max_oxo", 1)
    max_ooh_count = max_mod_cfg.get("max_ooh", 1)
    max_epoxy_count = max_mod_cfg.get("max_ep", 0)

    act_mod = sum([act_oh_count, act_oxo_count, act_epoxy_count, act_ooh_count])
    act_o = act_mod + act_ooh_count

    if all(
        [
            act_mod <= max_mod_count,
            act_o <= max_o_count,
            act_oh_count <= max_oh_count,
            act_oxo_count <= max_oxo_count,
            act_ooh_count <= max_ooh_count,
            act_epoxy_count <= max_epoxy_count,
            act_cho_count <= 1,
            act_cooh_count <= 1,
        ]
    ):
        is_allowed = True
    else:
        # print(f"Modification not possible: {smi}")
        pass

    return is_allowed


def predict_oap_segments(oap_smiles: list, unit_count: int = 1):
    mod_comb_tp_lst = list(itertools.product(oap_smiles, repeat=unit_count))
    # oap_mod_lst = ["".join(mod_units) for mod_units in mod_comb_tp_lst]  # all comb with position
    oap_mod_lst = natsorted(set(["".join(natsorted(list(m))) for m in mod_comb_tp_lst]))
    return oap_mod_lst


def predict_ocp_segments(ocp_smiles: list, oap_smiles: list, unit_count: int = 1):
    ocp_comb_lst = ocp_smiles.copy()
    if unit_count > 1:
        tmp_unit_count_lst = range(unit_count - 1, 0, -1)
        for tmp_unit_count in tmp_unit_count_lst:
            mod_comb_tp_lst = list(itertools.product(oap_smiles, repeat=tmp_unit_count))
            # mod_comb_lst = ["".join(mod_units) for mod_units in mod_comb_tp_lst]  # all comb with position
            tmp_pre_oap_comb_lst = natsorted(
                set(["".join(natsorted(list(m))) for m in mod_comb_tp_lst])
            )
            tmp_ocp_comb_lst = []
            for ocp_smi in ocp_smiles:
                tmp_pre_ocp_comb_lst = [
                    f"{seg_oap}{ocp_smi}" for seg_oap in tmp_pre_oap_comb_lst
                ]
                tmp_ocp_comb_lst.extend(tmp_pre_ocp_comb_lst)
            ocp_comb_lst.extend(tmp_ocp_comb_lst)
    else:
        pass

    return ocp_comb_lst


def assign_abbr(
    smi: str,
    max_mod_count: int = 2,
    max_o_count: int = 2,
    max_oh_count: int = 2,
    max_keto_count: int = 1,
    max_ep_count: int = 1,
    max_ooh_count: int = 1,
    ocp: bool = False,
):

    if re.match(r"^(?P<head>OC[(])([\\/\=CO\(\)\d])+(?P<tail>[)]=O)$", smi):
        res_typ = "FA"
    elif re.match(r"(OC)([\\/\=C])+C", smi):
        res_typ = "O-"
    elif re.match(r"(O\\C=C/)([\\/\=C])+C", smi):
        res_typ = "P-"
    else:
        res_typ = ""

    c_count = smi.count("C")
    db_count = len(re.findall(r"C=C", smi))

    res_base_str = f"{res_typ}{c_count}:{db_count}"

    oh_count = len(re.findall(r"C\(O\)", smi))
    oxo_count = len(re.findall(r"C\(=O\)", smi))
    ooh_count = len(re.findall(r"C\(OO\)", smi))
    ep_count = len(re.findall(r"C1C\(O1\)C", smi))

    mod_count = oh_count + oxo_count + ooh_count + ep_count
    o_count = oh_count + oxo_count + 2 * ooh_count + ep_count

    if re.match(r".*C=O\)=O$", smi):
        oxo_count += 1

    if re.match(r".*C\(O\)=O\)=O$", smi):
        oh_count -= 1
        if oh_count < 0:
            oh_count = 0

    mod_order_lst = [
        (oh_count, "OH"),
        (oxo_count, "oxo"),
        (ooh_count, "OOH"),
        (ep_count, "Ep"),
    ]
    mod_seg_lst = []
    for mod_seg in mod_order_lst:
        if mod_seg[0] > 1:
            mod_seg_lst.append(f"{mod_seg[0]}{mod_seg[1]}")
        elif mod_seg[0] == 1:
            mod_seg_lst.append(f"{mod_seg[1]}")
        else:
            pass

    if re.match(r".*C\(O\)=O\)=O$", smi):
        mod_seg_lst.append("COOH")

    if mod_seg_lst:
        mod_seg_str = "<" + ",".join(mod_seg_lst) + ">"
    else:
        mod_seg_str = ""
    res_str = res_base_str + mod_seg_str

    if ocp and "COOH" not in mod_seg_lst:
        max_keto_count += 1

    if (
        mod_count <= max_mod_count
        and o_count <= max_o_count
        and oh_count <= max_oh_count
        and oxo_count <= max_keto_count
        and ep_count <= max_ep_count
        and ooh_count <= max_ooh_count
    ):
        pass
    else:
        res_str = None

    return res_str


def predict_one(
    mod_info: pd.DataFrame,
    smi: str = r"OC(CCCCCCC\C=C/C\C=C/CCCCC)=O",
    max_mod_count: int = 2,
    max_o_count: int = 2,
    max_oh_count: int = 2,
    max_keto_count: int = 1,
    max_ep_count: int = 1,
    max_ooh_count: int = 1,
):
    # res_mod_dct = {}
    res_mod_lite_dct = {}
    oap_seg_smi_lst = natsorted(
        mod_info[mod_info["OAP"] == 1]["SMILES"].values.tolist()
    )
    ocp_seg_smi_lst = natsorted(
        mod_info[mod_info["OCP"] == 1]["SMILES"].values.tolist()
    )
    # smi_df = mod_info.set_index("SMILES")
    # smi_dct = smi_df.to_dict(orient="index")

    max_mod_dct = {
        "max_mod": max_mod_count,
        "max_o": max_o_count,
        "max_oh": max_oh_count,
        "max_oxo": max_keto_count,
        "max_ooh": max_ooh_count,
        "max_ep": max_ep_count,
    }

    is_op_link = False
    if re.match(r"(OC|O\\C=C/)([\\/\=C\d])+C", smi):
        is_op_link = True

    re_match = res_rgx.match(smi)
    if re_match and not is_op_link:
        m_groups_dct = re_match.groupdict()
        # {'head': 'OC(', 'left_c': 'CCCCCC', 'units': 'C\\C=C/C\\C=C/', 'right_c': 'CCCCC', 'tail': ')=O'}
        res_head = m_groups_dct.get("head", r"OC(")
        left_c = m_groups_dct.get("left_c", "")
        units = m_groups_dct.get("units", r"C\C=C//")
        right_c = m_groups_dct.get("right_c", "")
        tail = m_groups_dct.get("tail", r")=O")

        unit_count = len(unit_rgx.findall(units))
        if unit_count > 0:
            oap_comb_lst = predict_oap_segments(oap_seg_smi_lst, unit_count)
            ocp_comb_lst = predict_ocp_segments(
                ocp_seg_smi_lst, oap_seg_smi_lst, unit_count
            )
            chk_oap_comb_lst = [
                f"{res_head}{left_c}{a}{right_c}{tail}"
                for a in oap_comb_lst
                if fit_mod_count(a, max_mod_dct)
            ]
            chk_ocp_comb_lst = [
                f"{res_head}{left_c}{c}{tail}" for c in ocp_comb_lst if fit_mod_count(c)
            ]

            oap_lite_dct = {}
            for res_mod_smi in chk_oap_comb_lst:
                res_mod_abbr = assign_abbr(
                    res_mod_smi,
                    max_mod_count,
                    max_o_count,
                    max_oh_count,
                    max_keto_count,
                    max_ep_count,
                    max_ooh_count,
                )
                if res_mod_abbr is not None and res_mod_abbr not in oap_lite_dct:
                    oap_lite_dct[res_mod_abbr] = {
                        "name": res_mod_abbr,
                        # "lipid_class": "residue",
                        "mod_type": "oap",
                        "smiles": res_mod_smi,
                    }
                    oap_lite_dct[res_mod_abbr].update(
                        calc_smi(lipid_class="SN", smiles=res_mod_smi)
                    )
            ocp_lite_dct = {}
            for res_mod_smi in chk_ocp_comb_lst:
                res_mod_abbr = assign_abbr(
                    res_mod_smi, max_mod_count, max_o_count, max_keto_count, ocp=True
                )
                if res_mod_abbr is not None and res_mod_abbr not in ocp_lite_dct:
                    ocp_lite_dct[res_mod_abbr] = {
                        "name": res_mod_abbr,
                        # "lipid_class": "residue",
                        "mod_type": "ocp",
                        "smiles": res_mod_smi,
                    }
                    ocp_lite_dct[res_mod_abbr].update(
                        calc_smi(lipid_class="SN", smiles=res_mod_smi)
                    )
            res_mod_lite_dct.update(oap_lite_dct)
            res_mod_lite_dct.update(ocp_lite_dct)

        else:
            print(f"SMILES not used for oxidation: {smi}")
    else:
        print(f"SMILES not used for oxidation: {smi}")

    return res_mod_lite_dct


def predict_residues(
    mod_info: pd.DataFrame,
    res_info: dict,
    max_mod_count: int = 2,
    max_o_count: int = 2,
    max_oh_count: int = 2,
    max_keto_count: int = 1,
    max_ep_count: int = 1,
    max_ooh_count: int = 1,
):

    mod_res_dct = {}
    for res in res_info:
        res_dct = res_info[res]
        smi_key = None
        for smi_k in ["smiles", "smiles", "SMILES"]:
            if smi_k in res_dct:
                smi_key = smi_k
        if smi_key:
            res_smi = res_info[res].get(smi_key)
            mod_res = predict_one(
                mod_info,
                smi=res_smi,
                max_mod_count=max_mod_count,
                max_o_count=max_o_count,
                max_oh_count=max_oh_count,
                max_keto_count=max_keto_count,
                max_ep_count=max_ep_count,
                max_ooh_count=max_ooh_count,
            )
            mod_res_dct[res] = {
                "pr_name": res,
                "pr_smi": res_smi,
                "predictions": mod_res,
            }
        else:
            print("Can NOT find SMILES information.")

    return mod_res_dct


def sum_mod_res(res_mod_info: dict):
    predictions_dct = {}
    for abbr in res_mod_info:
        abbr_info_dct = res_mod_info[abbr]
        pr_abbr = abbr_info_dct.get("pr_name")
        abbr_product_dct = abbr_info_dct.get("predictions", {})
        if abbr_product_dct:
            for product_abbr in abbr_product_dct:
                if product_abbr not in predictions_dct:
                    abbr_product_dct[product_abbr].update({"origins": [pr_abbr]})
                    predictions_dct[product_abbr] = abbr_product_dct[product_abbr]
                else:
                    origins = predictions_dct[product_abbr].get("origins", [])
                    origins.append(pr_abbr)
                    predictions_dct[product_abbr]["origins"] = origins

    return predictions_dct


def get_epilipid_smi(lipid_class, fa_dct):

    if lipid_class in ["LPA", "LPC", "LPE", "LPG", "LPI", "LPS"]:
        if lipid_class.upper() in list(pl_hg_dct.keys()):
            pl_hg = pl_hg_dct[lipid_class.upper()]
            gly_part = r")C"
            # fa1 FA 16:0
            # fa1 = r'OC(CCCCCCCCCCCCCCC)=O'
            # fa2 FA 18:2 (9Z, 12Z)
            # fa2 = r'OC(CCCCCCC/C=C\C/C=C\CCCCC)=O'
            pl_end = r")=O"
            # the following order is important!!
            smi_str = "".join(
                [
                    pl_hg,
                    "O",
                    gly_part,
                    fa_dct["FA1"].get("smiles", "O"),
                    pl_end,
                ]
            )
        else:
            smi_str = ""
    elif lipid_class in ["CE"]:
        ce_left_smi = r"[C@]12(CC=C3C[C@@H]("
        ce_right_smi = r")CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@@]([H])([C@@](C)([H])CCCC(C)C)CC[C@@]21[H])[H]"

        smi_str = "".join(
            [
                ce_left_smi,
                fa_dct["FA1"].get("smiles", "O"),
                ce_right_smi,
            ]
        )

    elif lipid_class in ["PL", "PA", "PC", "PE", "PG", "PI", "PS"]:
        if lipid_class.upper() in list(pl_hg_dct.keys()):
            pl_hg = pl_hg_dct[lipid_class.upper()]
            gly_part = r")C"
            # fa1 FA 16:0
            # fa1 = r'OC(CCCCCCCCCCCCCCC)=O'
            # fa2 FA 18:2 (9Z, 12Z)
            # fa2 = r'OC(CCCCCCC/C=C\C/C=C\CCCCC)=O'
            pl_end = r")=O"
            # the following order is important!!
            smi_str = "".join(
                [
                    pl_hg,
                    fa_dct["FA2"].get("smiles", "O"),
                    gly_part,
                    fa_dct["FA1"].get("smiles", "O"),
                    pl_end,
                ]
            )
        else:
            smi_str = ""

    elif lipid_class in ["TG"]:
        gly_part_a = r"C(C"
        gly_part_b = r")("
        gly_part_c = r")C"
        # the following order is important!!
        # e.g. '[H]C(COC(CCC)=O)(OC(CCCCCCCCCCCCC)=O)COC(CCCCC)=O'
        # {'FA1': r'OC(CCC)=O', 'FA2': r'OC(CCCCCCCCCCCCC)=O', 'FA3': r'OC(CCCCC)=O'}
        if "" not in [
            fa_dct["FA1"].get("smiles", "O"),
            fa_dct["FA2"].get("smiles", "O"),
            fa_dct["FA3"].get("smiles", "O"),
        ]:
            smi_str = "".join(
                [
                    gly_part_a,
                    fa_dct["FA1"].get("smiles", "O"),
                    gly_part_b,
                    fa_dct["FA2"].get("smiles", "O"),
                    gly_part_c,
                    fa_dct["FA3"].get("smiles", "O"),
                ]
            )
        else:
            smi_str = ""
    else:
        smi_str = ""

    return smi_str


def get_lipids_comb(
    lipid_classes: list,
    residues: dict,
    modifications: pd.DataFrame,
    max_mod_count: int = 2,
    max_o_count: int = 2,
    max_oh_count: int = 2,
    max_keto_count: int = 1,
    max_ep_count: int = 1,
    max_ooh_count: int = 1,
):
    predicted_res_mod_info = predict_residues(
        mod_info=modifications,
        res_info=residues,
        max_mod_count=max_mod_count,
        max_o_count=max_o_count,
        max_oh_count=max_oh_count,
        max_keto_count=max_keto_count,
        max_ep_count=max_ep_count,
        max_ooh_count=max_ooh_count,
    )
    predicted_mod_res_info = sum_mod_res(predicted_res_mod_info)

    sum_res_info = residues.copy()
    sum_res_info.update(predicted_mod_res_info)

    unmod_res_lst = list(residues.keys())
    mod_res_lst = list(predicted_mod_res_info.keys())
    sum_res_lst = unmod_res_lst + mod_res_lst

    mono_res_classes_lst = []
    dual_res_classes_lst = []
    tri_res_classes_lst = []

    for lc in lipid_classes:
        if lc in ["LPA", "LPC", "LPE", "LPG", "LPI", "LPS", "CE"]:
            mono_res_classes_lst.append(lc)
        elif lc in ["DG", "PA", "PC", "PE", "PG", "PI", "PS"]:
            dual_res_classes_lst.append(lc)
        elif lc in ["TG"]:
            tri_res_classes_lst.append(lc)
        else:
            pass
    mono_res_comb_lst = list(itertools.product(mono_res_classes_lst, sum_res_lst))
    dual_res_comb_lst = list(
        itertools.product(dual_res_classes_lst, unmod_res_lst, mod_res_lst)
    )
    tri_res_comb_lst = list(
        itertools.product(
            tri_res_classes_lst, unmod_res_lst, unmod_res_lst, mod_res_lst
        )
    )

    comb_seg_lst = mono_res_comb_lst + dual_res_comb_lst + tri_res_comb_lst

    return comb_seg_lst, sum_res_info


def get_lipid_product_info(job, sum_res_info, frag_patterns) -> dict:

    predicted_epilipidome = {}
    for seg_tp in job:
        tmp_lipid_class = seg_tp[0]
        tmp_lipid_res_seg_lst = []
        tmp_res_dct = {"lipid_class": tmp_lipid_class}
        mod_type = "unmod"
        mod_pr_lst = []
        for res_sn_i in range(1, len(seg_tp)):
            tmp_res_i_abbr = seg_tp[res_sn_i]
            tmp_res_i_dct = sum_res_info.get(tmp_res_i_abbr).copy()
            tmp_res_dct[f"FA{res_sn_i}"] = tmp_res_i_dct
            tmp_lipid_res_seg_lst.append(tmp_res_i_abbr.strip("FA"))
            if tmp_res_i_dct.get("mod_type") != "unmod":
                mod_type = tmp_res_i_dct.get("mod_type")
            mod_pr_lst.append(tmp_res_i_dct.get("origins"))
        tmp_lipid_name = f'{tmp_lipid_class}({"_".join(tmp_lipid_res_seg_lst)})'
        tmp_lipid_dct = {"name": tmp_lipid_name}
        tmp_lipid_smi = get_epilipid_smi(tmp_lipid_class, tmp_res_dct)
        # tmp_lipid_dct["lipid_class"] = tmp_lipid_class
        tmp_lipid_dct["mod_type"] = mod_type
        tmp_lipid_dct["smiles"] = tmp_lipid_smi
        tmp_lipid_dct.update(
            calc_smi(lipid_class=tmp_lipid_class, smiles=tmp_lipid_smi)
        )
        tmp_lipid_dct["info"] = tmp_res_dct
        mod_pr_lst = [mod_pr for mod_pr in mod_pr_lst if mod_pr is not None]
        tmp_pr_comb_lst = list(itertools.product(*mod_pr_lst))  # List[Tuple[str]]
        tmp_pr_lst = []
        for tmp_pr_comb in tmp_pr_comb_lst:
            tmp_pr_res_lst = [t_res.strip("FA") for t_res in tmp_pr_comb]
            tmp_pr_name = f'{tmp_lipid_class}({"_".join(tmp_pr_res_lst)})'
            tmp_pr_lst.append(tmp_pr_name)
        tmp_lipid_dct["origins"] = tmp_pr_lst
        tmp_lipid_dct = generate_fragments(tmp_lipid_dct, frag_patterns=frag_patterns)

        predicted_epilipidome[tmp_lipid_name] = tmp_lipid_dct

    return predicted_epilipidome


def predict_epilipidome(
    lipid_classes: list,
    residues: dict,
    modifications: pd.DataFrame,
    frag_patterns: dict,
    max_mod_count: int = 2,
    max_o_count: int = 2,
    max_oh_count: int = 2,
    max_keto_count: int = 1,
    max_ep_count: int = 1,
    max_ooh_count: int = 1,
    worker: int = 4,
):

    predicted_epilipidome = {}
    comb_seg_lst, sum_res_info = get_lipids_comb(
        lipid_classes,
        residues,
        modifications,
        max_mod_count=max_mod_count,
        max_o_count=max_o_count,
        max_oh_count=max_oh_count,
        max_keto_count=max_keto_count,
        max_ep_count=max_ep_count,
        max_ooh_count=max_ooh_count,
    )
    comb_count = len(comb_seg_lst)
    if comb_count > 1:
        pass
    else:
        comb_count = 1
    if isinstance(worker, int):
        if worker > comb_count:
            worker = comb_count
        else:
            worker = min(comb_count, worker)
            if worker > 16:
                worker = 16
    else:
        worker = min(comb_count, 16)
    if worker > 1:
        pool = Pool(processes=worker)

        job_steps = math.ceil(comb_count / worker)
        job_split_idx_lst = list(range(0, comb_count, job_steps))[:-1]
        # job_split_idx_lst.append(comb_count-1)
        print("split jobs at index", job_split_idx_lst)
        job_lst = []
        for idx in range(1, len(job_split_idx_lst)):
            job_lst.append(
                comb_seg_lst[job_split_idx_lst[idx - 1] : job_split_idx_lst[idx]]
            )
        job_lst.append(
            comb_seg_lst[job_split_idx_lst[-1] :]
        )  # add the last job segment
        print(
            f"Split #{comb_count} tasks to {worker} workers each has {job_steps} tasks"
        )
        print("start multiprocessing...")
        result_lst = []
        for job in job_lst:
            r = pool.apply_async(
                get_lipid_product_info, args=(job, sum_res_info, frag_patterns)
            )  # dict
            result_lst.append(r)
        pool.close()
        pool.join()
        print("multiprocessing finished...")
        print("start merging results...")
        for r_dct in result_lst:
            if isinstance(r_dct, dict):
                result_dct = r_dct
            else:
                result_dct = r_dct.get()  # debug only
                # try:
                #     result_dct = r_dct.get()
                # except (KeyError, SystemError, ValueError):
                #     print(f"Cannot load results from multiprocessing...")
                #     result_dct = {}
            predicted_epilipidome.update(result_dct)
        print("results merged...")
    else:
        print("start processing...")
        predicted_epilipidome = get_lipid_product_info(
            comb_seg_lst, sum_res_info, frag_patterns
        )
        print("processing finished...")

    return predicted_epilipidome


def get_predictions(params):
    lipid_class_lst = params.get("lipid_class_lst")
    fa_cfg = params.get("input_fa_path")
    refine_epilipidome_path = params.get("refine_epilipidome_path", None)
    output_epilipidome_path = params.get("output_epilipidome_path")
    max_site = params.get("max_site")
    max_o = params.get("max_o")
    max_oh = params.get("max_oh")
    max_keto = params.get("max_keto")
    max_ooh = params.get("max_ooh")
    max_ep = params.get("max_ep")
    mod_cfg = params.get("mod_cfg")
    frag_cfg = params.get("frag_cfg")
    worker = params.get("worker_count", 2)

    # fa_cfg = r"../config/1_FA_list.csv"
    # mod_cfg = r"../config/2_DB_Mod_cfg_lv1.csv"
    # frag_cfg = r"../config/4_FragmentationPatterns.xlsx"

    res_dct = load_residues(fa_cfg)
    mod_info = load_modifications(mod_cfg)
    frag_patterns = load_frag_pattern(frag_cfg)

    usr_epilipidome = predict_epilipidome(
        lipid_classes=lipid_class_lst,
        residues=res_dct,
        modifications=mod_info,
        frag_patterns=frag_patterns,
        max_mod_count=max_site,
        max_o_count=max_o,
        max_oh_count=max_oh,
        max_keto_count=max_keto,
        max_ep_count=max_ep,
        max_ooh_count=max_ooh,
        worker=worker,
    )
    epilipidome_size = len(list(usr_epilipidome.keys()))
    epilipidome_size_info_str = f"Predicted epilipidome size:{epilipidome_size}"
    if isinstance(refine_epilipidome_path, str) and os.path.isfile(refine_epilipidome_path):
        usr_epilipidome = refine_epilipidome(usr_epilipidome, refine_epilipidome_path)
        epilipidome_size = len(list(usr_epilipidome.keys()))
        epilipidome_size_info_str += f"\nRefined epilipidome size:{epilipidome_size}"
    print(epilipidome_size_info_str)

    t1 = time.time()
    print("start to export file ...")
    with open(output_epilipidome_path, mode="w", encoding="utf-8") as fp:
        json.dump(usr_epilipidome, fp)
    print("file created...")

    run_info = "\n".join(
        [
            # f'Prediction finished in {time.time() - t1} sec',
            epilipidome_size_info_str,
            f"Output file saved.",
        ]
    )

    return run_info


# def get_multi_predictions(params: dict):
#
#     lipid_output_lst = []
#     fa_output_lst = []
#
#     for lipid_class in params:
#         info_updater_1, info_updater_2 = run_prediction(params[lipid_class])
#         print(info_updater_1)
#         print(info_updater_2)
#         output_path = params[lipid_class]["sdf_path"][:-4]
#         lipid_output_path = f'{output_path}.xlsx'
#         fa_output_path = f'{output_path}_FA_SUM.xlsx'
#         lipid_output_lst.append(lipid_output_path)
#         fa_output_lst.append(fa_output_path)
#
#     if lipid_output_lst:
#         lipid_df = merge_tables(lipid_output_lst)
#     else:
#         lipid_df = pd.DataFrame()
#     if fa_output_lst:
#         fa_df = merge_tables(fa_output_lst)
#     else:
#         fa_df = pd.DataFrame()
#
#     sum_lipid_output_path = f'{temp_file_path}_Merged.xlsx'
#     sum_fa_output_path = f'{temp_file_path}_Merged_FA_SUM.xlsx'
#
#     lipid_df.to_excel(sum_lipid_output_path)
#     fa_df.to_excel(sum_fa_output_path)
#
#     return {"predicted_lipids": lipid_df, "predicted_residues": fa_df}
#
#
# def predict_residues(names: List[str], config: dict):
#     saturated_fa_dct, unsaturated_fa_dct = filter_residues(names)
#     print(config)
#     for fa_name in unsaturated_fa_dct:
#         print(unsaturated_fa_dct[fa_name])


if __name__ == "__main__":

    max_worker = multiprocessing.cpu_count() - 1

    usr_params = {
        # "lipid_class_lst": ["PC"],
        "lipid_class_lst": ["PC", "PE"],
        # "lipid_class_lst": ["PA", "PC", "PE", "PG", "PS", "TG", "CE"],
        # "lipid_class_lst": ["TG"],
        # "lipid_class_lst": ["CE"],
        # "input_fa_path": r"config/1_FA_list.csv",
        "input_fa_path": r"config/1_FA_list_lite.csv",
        # "input_fa_path": r"config/1_FA_list_lite2.csv",
        # "refine_epilipidome_path": r"test/input/unmodified_lipids_PLs.xlsx",
        "output_epilipidome_path": r"test/output/oxPCPE_Lite.json",
        # "output_epilipidome_path": r"test/output/oxTG_lite.json",
        # "output_epilipidome_path": r"test/output/oxPL_TG_CE_mega.json",
        # "max_site": 2,
        # "max_o": 3,
        # "max_oh": 2,
        # "max_keto": 1,
        # "max_ooh": 1,
        # "max_ep": 1,
        "max_site": 2,
        "max_o": 3,
        "max_oh": 2,
        "max_keto": 1,
        "max_ooh": 1,
        "max_ep": 0,
        "mod_cfg": r"config/2_DB_Mod_cfg_lv1.csv",
        "frag_cfg": r"config/4_FragmentationPatterns.xlsx",
        "worker_count": max_worker
,
    }

    usr_run_info = get_predictions(usr_params)
    print(usr_run_info)
    print("FIN")
