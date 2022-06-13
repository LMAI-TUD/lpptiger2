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
import os
from typing import Dict, List

import pandas as pd
import pymzml


def extract_mzml_spectra(
    mzml: str,
    rt_range: List[float] = None,
    dda_top: int = 6,
    ms1_threshold: int = 1000,
    ms2_threshold: int = 10,
    ms1_precision: float = 50e-6,
    ms2_precision: float = 500e-6,
    min_spec_peaks: int = 3,
    vendor: str = "thermo",
    ms1_max: int = 0,
) -> (pd.DataFrame, Dict[int, pd.DataFrame], pd.DataFrame):
    """
    Extract mzML to a scan info DataFrame and a pandas panel for reader DataFrame of mz and i
    pymzml 2.2.5 is used

    Args:
        mzml (str): the file path of mzML file
        rt_range (list): a List of RT in minutes. e.g. [15, 30] for 15 to 30 min
        dda_top (int): DDA settings e.g. DDA TOP 6
        ms1_threshold (int): the absolute threshold for MS1 reader
        ms2_threshold (int): an absolute threshold for MS2 reader
        ms1_precision (float): e.g. 50e-6 for 50 ppm
        ms2_precision (float): e.g. 500e-6 for 500 ppm
        min_spec_peaks (int): minimum peaks a spectrum must have to be used for hunter, default = 3
        vendor (str): MS vendor abbreviations use lower case in list ['agilent', 'bruker', 'sciex', 'thermo', 'waters']
        ms1_max (int): Max of MS1 intensity, used to search for low intensity signals, set 0 to disable by default

    Returns:
        scan_info_df (pd.DataFrame):
        spec_dct (Dict[int, pd.DataFrame]):
        ms1_xic_df (pd.DataFrame):

    """

    if isinstance(rt_range, list) and len(rt_range) == 2:
        pass
    else:
        rt_range = [-1, 999]

    rt_start = rt_range[0]
    rt_end = rt_range[1]
    if os.path.isfile(mzml):
        print("[STATUS] >>> Start to process file: %s" % mzml)
    else:
        print("[ERROR] !!! FileNotFoundError: %s" % mzml)
        raise FileNotFoundError
    print(
        "[INFO] --> Processing RT: %.2f -> %.2f with DDA Top % i"
        % (rt_start, rt_end, dda_top)
    )
    try:
        spec_obj = pymzml.run.Reader(
            mzml, MS1_Precision=ms1_precision, MSn_Precision=ms2_precision
        )
    except (IOError, OSError):
        try:
            # Try to use the CV from mzML 4.0 released in 2016
            spec_obj = pymzml.run.Reader(
                mzml,
                MS1_Precision=ms1_precision,
                MSn_Precision=ms2_precision,
                obo_version="4.0.1",
            )
        except (IOError, OSError):
            # try to use legacy version 1.1.0 to parse the mzML
            spec_obj = pymzml.run.Reader(
                mzml,
                MS1_Precision=ms1_precision,
                MSn_Precision=ms2_precision,
                obo_version="1.1.0",
            )

    spec_idx = 1
    dda_event_idx = 0
    spec_idx_lst = []
    dda_event_lst = []
    rt_lst = []
    dda_rank_lst = []
    scan_id_lst = []
    pr_mz_lst = []
    polarity_lst = []

    scan_info_dct = {
        "spec_index": spec_idx_lst,
        "scan_time": rt_lst,
        "dda_event_idx": dda_event_lst,
        "DDA_rank": dda_rank_lst,
        "scan_number": scan_id_lst,
        "MS2_PR_mz": pr_mz_lst,
        "polarity": polarity_lst,
    }

    spec_dct = {}  # Type: Dict[pd.DataFrame]
    if vendor == "waters":
        ms2_function_range_lst = list(range(2, dda_top + 1))
    else:
        ms2_function_range_lst = [2]

    ms1_xic_df = pd.DataFrame()

    print("Instrument vendor: %s" % vendor)

    if vendor in ["agilent", "bruker", "sciex", "thermo", "waters"]:
        dda_rank_idx = 0
        for _spectrum in spec_obj:  # type: pymzml.spec.Spectrum

            # # Reserved legacy code for Ion mobility capabilities.
            # if ims_obo in list(_spectrum.keys()):
            #     try:
            #         if len(_spectrum.peaks) > 4:
            #             not_empty_spec = 1
            #         else:
            #             not_empty_spec = 0
            #     except (KeyError, ValueError):
            #         not_empty_spec = 0
            # else:
            #     not_empty_spec = 1

            pr_mz = 0
            try:
                _scan_rt = float(_spectrum.scan_time[0])
                if isinstance(_spectrum.scan_time[1], str) and _spectrum.scan_time[
                    1
                ].lower() in ["s", "sec", "second", "seconds"]:
                    _scan_rt = round(_scan_rt / 60, 6)
            except (ValueError, TypeError):
                _scan_rt = -0.1

            if (
                rt_start <= _scan_rt <= rt_end
                and _spectrum.mz.any()
                and _spectrum.id_dict
            ):
                try:
                    _scan_id = int(_spectrum.id_dict.get("scan", -1))
                except ValueError:
                    _scan_id = -1
                if _scan_id == -1:
                    try:
                        _scan_id = int(_spectrum.id_dict.get("ID", spec_idx))
                    except ValueError:
                        _scan_id = spec_idx
                try:
                    ms_level = int(_spectrum.ms_level)
                except ValueError:
                    print(_spectrum.ms_level)
                    ms_level = -1

                _raw_tmp_spec_df = pd.DataFrame(
                    data={"mz": _spectrum.mz, "i": _spectrum.i}
                )
                _tmp_spec_df = pd.DataFrame()
                polarity = ""
                is_negative_polarity = _spectrum["negative scan"]
                if is_negative_polarity is None:
                    pass
                else:
                    if is_negative_polarity is True:
                        polarity = "neg"
                is_positive_polarity = _spectrum["positive scan"]
                if is_positive_polarity is None:
                    pass
                else:
                    if is_positive_polarity is True:
                        polarity = "pos"
                if ms_level == 1 and _scan_id > 0:
                    dda_event_idx += 1  # a new set of DDA start from this new MS1
                    dda_rank_idx = (
                        0  # set the DDA rank back to 0 for the survey MS1 scans
                    )
                    # use ms1_threshold * 0.1 to keep isotope patterns
                    if ms1_max > ms1_threshold:
                        _tmp_spec_df = _raw_tmp_spec_df.query(
                            "%f <= i <= %f" % ((ms1_threshold * 0.1), ms1_max)
                        ).copy()
                    else:
                        _tmp_spec_df = _raw_tmp_spec_df.query(
                            "%f <= i" % (ms1_threshold * 0.1)
                        ).copy()

                    if not _tmp_spec_df.empty:
                        _tmp_spec_df = _tmp_spec_df.sort_values(by="i", ascending=False)
                        _tmp_spec_df = _tmp_spec_df.reset_index(drop=True)
                        _tmp_spec_df.loc[:, "rt"] = _scan_rt
                        ms1_xic_df = pd.concat([ms1_xic_df, _tmp_spec_df], ignore_index=True)
                    else:
                        print("empty_MS1_spectrum --> index = ", spec_idx)

                    print(
                        f"MS1 -> index: {spec_idx} @ scan_time: {_scan_rt:.3f} | mode: {polarity} | DDA_events: {dda_event_idx}"
                    )

                elif ms_level in ms2_function_range_lst and _scan_id > 0:
                    dda_rank_idx += 1
                    try:
                        pr_mz = _spectrum.selected_precursors[0].get("mz", -1)
                    except (KeyError, AttributeError):
                        pr_mz = -1
                    if pr_mz > 0:
                        _tmp_spec_df = _raw_tmp_spec_df.query("i >= %f" % ms2_threshold)
                        # print(_tmp_spec_df)
                        if _tmp_spec_df.shape[0] > min_spec_peaks:
                            pass
                        else:
                            print("empty_MS2_spectrum --> index = ", spec_idx)

                        print(
                            f"MS2 -> index: {spec_idx} @ scan_time: {_scan_rt:.3f} | mode: {polarity} "
                            f"| DDA_events: {dda_event_idx} | RANK {dda_rank_idx} | PR_MZ: {pr_mz}"
                        )
                    else:
                        print(
                            "MS2 DDA RANK ERROR of rank: {rank}".format(
                                rank=dda_rank_idx
                            )
                        )
                else:
                    print(
                        f"[ERROR] Can not read the spectrum # {spec_idx} @ {_scan_rt:.3f} min - ms_level {ms_level}"
                    )

                if not _tmp_spec_df.empty:
                    spec_dct[spec_idx] = _tmp_spec_df
                    spec_idx_lst.append(spec_idx)
                    dda_event_lst.append(dda_event_idx)
                    rt_lst.append(_scan_rt)
                    dda_rank_lst.append(dda_rank_idx)
                    scan_id_lst.append(_scan_id)
                    pr_mz_lst.append(pr_mz)
                    polarity_lst.append(polarity)

                del _raw_tmp_spec_df
            else:  # rt not in defined range, skip.
                pass
                # print(f'[INFO] Skip spectrum # {spec_idx} @ {_scan_rt:.3f} min - ms_level {ms_level}')
            spec_idx += 1

    else:
        raise ValueError(
            f"LipidHunter do not support mzML from this vendor: {vendor}\n"
            f'Supported vendors: ["agilent", "bruker", "sciex", "thermo", "waters"]'
        )
    scan_info_df = pd.DataFrame(
        data=scan_info_dct,
        columns=[
            "dda_event_idx",
            "spec_index",
            "scan_time",
            "DDA_rank",
            "scan_number",
            "MS2_PR_mz",
            "polarity",
        ],
    )
    scan_info_df.sort_values(by="scan_time", inplace=True)
    scan_info_df = scan_info_df.round({"MS2_PR_mz": 6})
    int_col_lst = ["dda_event_idx", "spec_index", "DDA_rank", "scan_number"]
    scan_info_df[int_col_lst] = scan_info_df[int_col_lst].astype(int)
    print("=== ==> --> mzML extracted")

    return scan_info_df, spec_dct, ms1_xic_df


def find_ms1_pr(ms1_df: pd.DataFrame, mz_lib: float, ms1_ppm: int):
    ms1_df = ms1_df.query("i > 0")
    ms1_df = ms1_df.sort_values(by="i", ascending=False).reset_index(drop=True)
    ms1_delta = mz_lib * ms1_ppm * 0.000001
    ms1_pr_df = ms1_df.query(
        f"{mz_lib - ms1_delta} <= mz <= {mz_lib + ms1_delta}"
    ).copy()

    if not ms1_pr_df.empty:
        ms1_pr_df["mz"] = ms1_pr_df["mz"].round(6)
        ms1_pr_df = ms1_pr_df.sort_values(by="i", ascending=False)
        ms1_pr_se = ms1_pr_df.iloc[0]
        ms1_obs_pr_mz = ms1_pr_se["mz"]
        ms1_obs_pr_i = ms1_pr_se["i"]
        ms1_obs_pr_ppm = (1e6 * (ms1_obs_pr_mz - mz_lib) / mz_lib).round(2)
    else:
        ms1_obs_pr_mz = 0
        ms1_obs_pr_i = 0
        ms1_obs_pr_ppm = 999.9
        # print("[WARNING] !!! Precursor m/z in MS1 not in the list ...")

    return ms1_obs_pr_mz, ms1_obs_pr_i, ms1_obs_pr_ppm


if __name__ == "__main__":

    test_mzml = r"../../temp/QE_21_24_MF49.mzML"
    scan_info_df, spec_dct, usr_ms1_xic_df = extract_mzml_spectra(test_mzml)
