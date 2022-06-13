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
import json
import random

from PySide6 import QtCore, QtGui
import matplotlib

# matplotlib.use('Qt4Agg')
# matplotlib.use("agg")
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import pandas as pd


def plot_xic(
    xic_fig,
    xic_info: dict,
    xic_title_str: str,
    ms1_scan_time: float,
    ms2_scan_time: float,
):
    # print('start to plot XIC ...')
    xic_peaks_dct = xic_info.get("peaks", [])
    if xic_peaks_dct:
        xic_peaks_df = pd.DataFrame(data=xic_peaks_dct, columns=["rt", "i"])
    else:
        xic_peaks_df = pd.DataFrame()
    if not xic_peaks_df.empty:
        xic_fig.tick_params(axis="both", which="major", labelsize=10)
        xic_rt_min = xic_peaks_df["rt"].min()
        xic_rt_max = xic_peaks_df["rt"].max()
        xic_i_max = xic_peaks_df["i"].max()
        xic_rt_label_shift = (xic_rt_max - xic_rt_min) * 0.04
        xic_fig.plot(xic_peaks_df["rt"], xic_peaks_df["i"], alpha=0.7, color="grey")
        xic_fig.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
        _marker_line, _stem_lines, _base_line = xic_fig.stem(
            [ms1_scan_time], [xic_i_max], markerfmt=" "
        )
        plt.setp(_stem_lines, color=(0.0, 0.4, 1.0, 0.3), linewidth=3)
        _marker_line, _stem_lines, _base_line = xic_fig.stem(
            [ms2_scan_time], [xic_i_max], linefmt="--", markerfmt=" "
        )
        plt.setp(_stem_lines, color=(0.0, 0.65, 1.0, 0.95), linewidth=2)
        xic_fig.text(
            ms1_scan_time - xic_rt_label_shift,
            xic_i_max * 0.98,
            "MS",
            fontsize=8,
            color=(0.0, 0.4, 1.0, 1.0),
            weight="bold",
        )
        xic_fig.text(
            ms2_scan_time,
            xic_i_max * 0.98,
            "MS/MS",
            fontsize=8,
            color=(0.0, 0.65, 1.0, 1.0),
            weight="bold",
        )
        xic_fig.set_xlabel("Scan time (min)", fontsize=7, labelpad=-1)
        xic_fig.set_ylabel("Intensity", fontsize=7)
        xic_fig.set_xlim([xic_rt_min, xic_rt_max])
        xic_fig.set_ylim([0, xic_i_max * 1.1])
        xic_fig.set_title(xic_title_str, color=(0.0, 0.4, 1.0, 1.0), fontsize=8, y=0.98)


def plot_ms1_fig(ms1_fig, ms1_info: dict, ms1_title_str: str, ms1_obs_mz, ms1_obs_i):
    # print('start to plot MS1 ...')
    ms1_fig.tick_params(axis="both", which="major", labelsize=10)
    ms1_df = ms1_info.get("df", pd.DataFrame())

    ms1_obs_mz = ms1_info["obs_mz"]
    ms1_obs_i = ms1_info["obs_i"]

    if ms1_df.shape[0] > 700:  # switch plot mode
        ms1_fig.plot(ms1_df["mz"], ms1_df["i"], "grey", lw=0.6)
    else:
        marker_l, stem_l, base_l = ms1_fig.stem(
            ms1_df["mz"], ms1_df["i"], markerfmt=" "
        )
        plt.setp(stem_l, color="grey", lw=0.6)
        plt.setp(base_l, visible=False)
    _marker_l, _stem_l, _base_l = ms1_fig.stem(
        [ms1_obs_mz], ms1_df["i"].max(), markerfmt=" "
    )
    plt.setp(_stem_l, color=(0.0, 0.4, 1.0, 0.2), linewidth=4)
    plt.setp(_base_l, visible=False)
    _marker_l, _stem_l, _base_l = ms1_fig.stem([ms1_obs_mz], [ms1_obs_i], markerfmt="D")
    # plt.setp(_stem_lines, color=(0.4, 1.0, 0.8, 1.0))
    plt.setp(
        _marker_l, markerfacecolor=(0.3, 0.9, 1.0, 0.8), markersize=5, markeredgewidth=0
    )
    plt.setp(_stem_l, visible=False)
    plt.setp(_base_l, visible=False)
    ms1_fig.text(
        ms1_obs_mz,
        ms1_obs_i * 1.036,
        f"{ms1_obs_mz:.4f}",
        fontsize=7,
        color=(0.0, 0.4, 1.0, 1.0),
        rotation=60,
    )
    ms1_fig.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    ms1_fig.set_xlabel("m/z", fontsize=7, labelpad=-1)
    ms1_fig.set_ylabel("Intensity", fontsize=7)
    ms1_fig.set_ylim([0, max(ms1_df["i"].values.tolist()) * 1.36])

    # add annotation
    pre_ms_pkl_top_df = ms1_df.sort_values(by="i", ascending=False).head(30)
    # avoid overlay with pr_mz
    _ms_pkl_top_df = pre_ms_pkl_top_df[
        ~pre_ms_pkl_top_df.iloc[:, 0].between(
            ms1_obs_mz - 3, ms1_obs_mz + 3, inclusive="both"
        )
    ]
    _ms_pkl_top_peak_list = list(
        zip(_ms_pkl_top_df["mz"].values.tolist(), _ms_pkl_top_df["i"].values.tolist())
    )
    skipped_mz_lst = []

    for _ms_pkl_top_peak in _ms_pkl_top_peak_list:
        ms1_label_mz = _ms_pkl_top_peak[0]
        if ms1_label_mz not in skipped_mz_lst:  # avoid overlay with previous plotted labels
            _ms_pkl_top_peak_str = f"{ms1_label_mz:.4f}"
            _ms_pkl_top_peak_y = _ms_pkl_top_peak[1]
            ms1_fig.text(
                ms1_label_mz,
                _ms_pkl_top_peak_y,
                _ms_pkl_top_peak_str,
                fontsize=6,
                rotation=60,
            )
            _remove_peaks_df = _ms_pkl_top_df[_ms_pkl_top_df.iloc[:, 0].between(
                    ms1_label_mz - 2, ms1_label_mz + 2, inclusive="both")]
            skipped_mz_lst.extend(_remove_peaks_df["mz"].to_list())

    ms1_fig.set_title(ms1_title_str, color=(0.0, 0.4, 1.0, 1.0), fontsize=8, y=0.98)
    # print(core_count, 'plot MS in ', time.time() - _t_img_0)


def plot_ms1_zoomed_fig(
    ms1_zoomed_fig,
    ms1_info: dict,
    ms_zoomed_title_str: str,
):
    # print('start to plot MS zoom ...')
    ms1_zoomed_fig.tick_params(axis="both", which="major", labelsize=10)
    ms1_df = ms1_info.get("df", pd.DataFrame())

    lib_mz = ms1_info["lib_mz"]
    ms1_obs_mz = ms1_info["obs_mz"]
    ms1_obs_i = ms1_info["obs_i"]
    ms1_delta = lib_mz * ms1_info.get("ms1_tolerance_ppm", 50) * 0.000001

    ms1_zoomed_df = ms1_df.query(f"{ms1_obs_mz} - 1.5 <= mz <= {ms1_obs_mz} + 3.55")

    ms1_zoomed_fig.set_xlim([ms1_obs_mz - 1.5, ms1_obs_mz + 3.55])
    ms_zoomed_bp_i = ms1_zoomed_df["i"].max()
    ms1_zoomed_fig.set_ylim([0, ms_zoomed_bp_i * 1.45])
    ms1_zoomed_fig.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    ms1_zoomed_fig.ticklabel_format(axis="x", useOffset=False)
    ms1_zoomed_fig.set_xlabel("m/z", fontsize=7, labelpad=-1)
    ms1_zoomed_fig.set_ylabel("Intensity", fontsize=7)

    plt_ms_zoom_df = ms1_zoomed_df.sort_values(by="mz", ascending=True)
    theo_i_bar_color = (0, 0.4, 1.0, 0.5)
    theo_i2_bar_color = (0, 0.4, 1.0, 0.3)

    # isotope region | if any peak in M-1.0034

    isotope_labels_dct = ms1_info.get("isotope_labels", {})
    isotope_signal_labels_dct = isotope_labels_dct.get("signal_labels", {})
    isotope_noise_labels_dct = isotope_labels_dct.get("noise_labels", {})
    ms_zoom_offset_i = ms_zoomed_bp_i * 0.15

    pseudo_m0_mz = 0
    pseudo_m0_i = 0
    if "m°" in isotope_noise_labels_dct:
        pseudo_m0_mz = isotope_noise_labels_dct.get("m°", [0, 0])[0]
        pseudo_m0_i = isotope_noise_labels_dct.get("m°", [0, 0])[1]

    if pseudo_m0_i:
        pseudo_m0_theo_i_box = patches.Rectangle(
            (lib_mz - 1.0034 - ms1_delta, 0),
            2 * ms1_delta,
            pseudo_m0_i,
            facecolor=(1.0, 0.0, 0.0, 0.5),
            edgecolor="none",
            zorder=9,
        )
        ms1_zoomed_fig.add_patch(pseudo_m0_theo_i_box)
        if pseudo_m0_i > 0.05 * ms1_obs_i:
            ms1_zoomed_fig.text(
                pseudo_m0_mz + 0.06,
                pseudo_m0_i,
                f"{pseudo_m0_mz:.4f}",
                color="red",
                fontsize=7,
                alpha=0.8,
            )
            ms1_zoomed_fig.text(
                pseudo_m0_mz - 0.06,
                pseudo_m0_i + ms_zoom_offset_i,
                "[m°]",
                color="red",
                fontsize=7,
                alpha=0.8,
            )
    m0_theo_base_i = isotope_noise_labels_dct.get("m°+1", [0, 0])[1]
    if m0_theo_base_i:
        m0_theo_base_box = patches.Rectangle(
            (lib_mz - ms1_delta, 0),
            2 * ms1_delta,
            m0_theo_base_i,
            facecolor=(1.0, 0.0, 0.0, 0.5),
            edgecolor="none",
            zorder=1,
        )
        ms1_zoomed_fig.add_patch(m0_theo_base_box)
    m0_theo_box = patches.Rectangle(
        (lib_mz - ms1_delta, m0_theo_base_i),
        2 * ms1_delta,
        ms1_obs_i - m0_theo_base_i,
        facecolor=theo_i_bar_color,
        edgecolor="none",
        zorder=10,
    )
    ms1_zoomed_fig.add_patch(m0_theo_box)

    _marker_l, _stem_l, _base_l = ms1_zoomed_fig.stem(
        [ms1_obs_mz], [ms1_obs_i], markerfmt="D"
    )
    plt.setp(_stem_l, visible=False)
    plt.setp(_base_l, visible=False)
    plt.setp(
        _marker_l,
        markerfacecolor=(0.2, 0.8, 1.0, 0.8),
        markeredgecolor="none",
        markeredgewidth=0,
        markersize=6,
        lw=1,
        zorder=21,
    )

    ms1_zoomed_fig.text(
        ms1_obs_mz + 0.06,
        ms1_obs_i,
        f"{ms1_obs_mz:.4f}",
        color=(0.0, 0.4, 1.0, 1.0),
        fontsize=7,
    )

    ms1_zoomed_fig.text(
        lib_mz - 0.2,
        ms1_obs_i + ms_zoom_offset_i,
        "[M]",
        color=(0.2, 0.8, 1.0, 1.0),
        fontsize=7,
    )
    ms1_zoomed_fig.text(
        lib_mz - 0.78,
        ms1_obs_i,
        f"Calc: {lib_mz:.4f}",
        color=(0.2, 0.8, 1.0, 1.0),
        fontsize=7,
    )

    # isotope region | highlight the 1st isotope
    # theo range box
    m1_theo_mz = isotope_signal_labels_dct.get("M+1", [0, 0])[0]
    m1_theo_i = isotope_signal_labels_dct.get("M+1", [0, 0])[1]
    m1_theo_base_i = isotope_noise_labels_dct.get("m°+2", [0, 0])[1]
    if m1_theo_base_i:
        m1_theo_base_box = patches.Rectangle(
            (m1_theo_mz - ms1_delta, 0),
            2 * ms1_delta,
            m1_theo_base_i,
            facecolor=(1.0, 0.0, 0.0, 0.5),
            edgecolor="none",
            zorder=11,
        )
        ms1_zoomed_fig.add_patch(m1_theo_base_box)

    m1_theo_box = patches.Rectangle(
        (m1_theo_mz - ms1_delta, m1_theo_base_i),
        2 * ms1_delta,
        m1_theo_i - m1_theo_base_i,
        facecolor=theo_i_bar_color,
        edgecolor="none",
        zorder=12,
    )
    ms1_zoomed_fig.add_patch(m1_theo_box)

    _marker_l, _stem_l, _base_l = ms1_zoomed_fig.stem(
        [m1_theo_mz], [m1_theo_i], linefmt="--", markerfmt="o"
    )
    plt.setp(
        _marker_l,
        markerfacecolor=(0.2, 0.8, 1.0, 0.8),
        markersize=6,
        markeredgewidth=0,
        zorder=22,
    )
    plt.setp(_stem_l, visible=False)
    plt.setp(_base_l, visible=False)

    m1_obs_df = ms1_zoomed_df.query(
        f"{m1_theo_mz - 2 * ms1_delta} <= mz <= {m1_theo_mz + 2 * ms1_delta}"
    ).copy()  # pd.DataFrame

    if not m1_obs_df.empty:
        m1_obs_df.sort_values(by="i", ascending=False)
        m1_obs_mz = m1_obs_df.iloc[0][0]
        ms1_zoomed_fig.text(
            m1_obs_mz + 0.06,
            m1_obs_df.iloc[0][1],
            f"{m1_obs_mz:.4f}",
            color=(0.0, 0.4, 1.0, 1.0),
            fontsize=7,
        )

    ms1_zoomed_fig.text(
        m1_theo_mz - 0.35,
        m1_theo_i + ms_zoom_offset_i,
        "[M+1]",
        color=(0.2, 0.8, 1.0, 1.0),
        fontsize=7,
    )

    ms1_zoomed_fig.text(
        m1_theo_mz - 0.78,
        m1_theo_i,
        f"Calc: {m1_theo_mz:.4f}",
        color=(0.2, 0.8, 1.0, 1.0),
        fontsize=7,
    )

    # # isotope region | highlight the 2nd isotope
    m2_theo_mz = isotope_signal_labels_dct.get("M+2", [0, 0])[0]
    m2_theo_i = isotope_signal_labels_dct.get("M+2", [0, 0])[1]
    # m2_theo_base_i = isotope_noise_labels_dct.get("m°+3", [0, 0])[1]

    if m2_theo_mz and m2_theo_i:

        m2_theo_box = patches.Rectangle(
            (m2_theo_mz - ms1_delta, 0),
            2 * ms1_delta,
            m2_theo_i,
            facecolor=theo_i_bar_color,
            edgecolor="none",
            zorder=13,
        )
        ms1_zoomed_fig.add_patch(m2_theo_box)
        _marker_l, _stem_l, _base_l = ms1_zoomed_fig.stem(
            [m2_theo_mz], [m2_theo_i], linefmt="--", markerfmt="o"
        )
        plt.setp(
            _marker_l,
            markerfacecolor=(0.2, 0.8, 1.0, 0.8),
            markersize=6,
            markeredgewidth=0,
            zorder=23,
        )
        plt.setp(_stem_l, visible=False)
        plt.setp(_base_l, visible=False)

        m2_obs_df = ms1_zoomed_df.query(
            f"{m2_theo_mz - 2 * ms1_delta} <= mz <= {m2_theo_mz + 2 * ms1_delta}"
        ).copy()  # pd.DataFrame

        if not m2_obs_df.empty:
            m2_obs_df.sort_values(by="i", ascending=False)
            m2_obs_mz = m2_obs_df.iloc[0][0]
            ms1_zoomed_fig.text(
                m2_obs_mz + 0.06,
                m2_obs_df.iloc[0][1],
                f"{m2_obs_mz:.4f}",
                color=(0.0, 0.4, 1.0, 1.0),
                fontsize=7,
            )

        ms1_zoomed_fig.text(
            m2_theo_mz - 0.35,
            m2_theo_i + ms_zoom_offset_i,
            "[M+2]",
            color=(0.2, 0.8, 1.0, 1.0),
            fontsize=7,
        )

        ms1_zoomed_fig.text(
            m2_theo_mz - 0.78,
            m2_theo_i,
            f"Calc: {m2_theo_mz:.4f}",
            color=(0.2, 0.8, 1.0, 1.0),
            fontsize=7,
        )

        # plot M+2
        mh2_theo_mz = isotope_noise_labels_dct.get("m'", [0, 0])[0]
        mh2_theo_i = isotope_noise_labels_dct.get("m'", [0, 0])[1]

        # mh2_theo_base_box = patches.Rectangle(
        #     (mh2_theo_mz - ms1_delta, 0),
        #     2 * ms1_delta,
        #     deconv_lst[decon_idx],
        #     facecolor=theo_i2_bar_color,
        #     edgecolor="none",
        #     zorder=14,
        # )
        # ms1_zoomed_fig.add_patch(mh2_theo_base_box)
        if mh2_theo_i:
            ms1_zoomed_fig.text(
                mh2_theo_mz - 0.7,
                mh2_theo_i + m2_theo_i + ms_zoom_offset_i,
                "[m'] predicted",
                color="pink",
                fontsize=6,
            )

            ms1_zoomed_fig.text(
                mh2_theo_mz - 0.65,
                mh2_theo_i + m2_theo_i,
                f"Calc: {mh2_theo_mz:.4f}",
                color="pink",
                fontsize=6,
            )

            mh2_theo_box = patches.Rectangle(
                (mh2_theo_mz - ms1_delta, m2_theo_i),
                2 * ms1_delta,
                mh2_theo_i,
                facecolor=(1.0, 0.0, 0.0, 0.5),
                edgecolor="none",
                zorder=15,
            )
            ms1_zoomed_fig.add_patch(mh2_theo_box)
            _marker_l, _stem_l, _base_l = ms1_zoomed_fig.stem(
                [mh2_theo_mz], [mh2_theo_i + m2_theo_i], linefmt="--", markerfmt="o"
            )
            plt.setp(_stem_l, color="pink", alpha=0.6)
            plt.setp(
                _marker_l,
                markerfacecolor="pink",
                markersize=5,
                markeredgewidth=0,
                alpha=0.6,
                zorder=24,
            )
            plt.setp(_stem_l, visible=False)
        # plot M+2 + 1
        mh21_theo_mz = isotope_noise_labels_dct.get("m'+1", [0, 0])[0]
        mh21_theo_i = isotope_noise_labels_dct.get("m'+1", [0, 0])[1]
        if mh21_theo_i:
            ms1_zoomed_fig.text(
                mh21_theo_mz - 0.3,
                mh21_theo_i + ms_zoom_offset_i,
                "[m'+1]",
                color="pink",
                fontsize=6,
            )
            ms1_zoomed_fig.text(
                mh21_theo_mz - 0.65,
                mh21_theo_i,
                f"Calc: {mh21_theo_mz:.4f}",
                color="pink",
                fontsize=6,
            )

        # plot the M+2H isotope score
        mh2_isotope_score = ms1_info.get("m'_isotope_score", 0)
        ms1_zoomed_fig.text(
            lib_mz + 2.85,
            max(ms_zoomed_bp_i - ms_zoom_offset_i, ms1_obs_i * 0.8),
            f"[m'] isotope score = {mh2_isotope_score:.1f}",
            verticalalignment="top",
            horizontalalignment="right",
            color="pink",
            fontsize=6,
        )

    # plot the isotope score
    ms1_zoomed_fig.text(
        lib_mz + 0.2,
        ms_zoomed_bp_i + 2.5 * ms_zoom_offset_i,
        f"isotope score = {ms1_info.get('isotope_score', 0):.1f}",
        verticalalignment="top",
        horizontalalignment="right",
        color=(0.0, 0.4, 1.0, 1.0),
        fontsize=10,
    )

    # plot ms_zoom peaks
    if plt_ms_zoom_df.shape[0] > 999:
        ms1_zoomed_fig.plot(
            plt_ms_zoom_df["mz"].values.tolist(),
            plt_ms_zoom_df["i"].values.tolist(),
            "grey",
            lw=1,
            zorder=1,
        )
    else:
        marker_l, stem_l, base_l = ms1_zoomed_fig.stem(
            plt_ms_zoom_df["mz"].values.tolist(),
            plt_ms_zoom_df["i"].values.tolist(),
            markerfmt=" ",
        )
        plt.setp(stem_l, color="grey", lw=0.75, alpha=0.6)
        plt.setp(stem_l, zorder=1)
        plt.setp(base_l, visible=False)

    ms1_zoomed_fig.set_title(
        ms_zoomed_title_str, color=(0.0, 0.4, 1.0, 1.0), fontsize=8, y=0.98
    )


def plot_ms2_fig(ms2_fig, ms2_info: dict, ms2_title_str: str):

    ms2_df = ms2_info.get("df", pd.DataFrame())
    base_peak_i = ms2_df["i"].max()
    # msp_factor_i = -0.001 * base_peak_i

    ms2_fig.tick_params(axis="both", which="major", labelsize=10)
    ms2_fig.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    ms2_fig.set_xlabel("m/z", fontsize=7, labelpad=-1)
    ms2_fig.set_ylabel("Intensity", fontsize=7)

    # plot hunter table
    ident_col_labels = ("Proposed_structure", "Score")
    _ident_table_df = pd.DataFrame(
        data={
            "Proposed_structure": ms2_info.get("lipid_name", ""),
            "Score": f"{ms2_info.get('scores', {}).get('total_score', 0):.1f}",
        },
        index=[1],
    )
    ident_table_vals = list(map(list, _ident_table_df.values))
    ident_col_width_lst = [0.6, 0.15]
    ident_table = ms2_fig.table(
        cellText=ident_table_vals,
        colWidths=ident_col_width_lst,
        colLabels=ident_col_labels,
        loc="upper center",
        cellLoc="center",
    )
    ident_table.set_fontsize(8)

    # plot rank marks
    # rank_labels_dct = ms2_info.get("rank_labels", {})

    # plot msp marks
    msp_signal_labels_dct = ms2_info.get("msp_labels", {}).get("signal_labels", {})
    if msp_signal_labels_dct:
        msp_signal_mz_lst = [
            msp_signal_labels_dct[msp_key][0] for msp_key in msp_signal_labels_dct
        ]
        msp_signal_i_lst = [
            msp_signal_labels_dct[msp_key][1] * -1
            # msp_signal_labels_dct[msp_key][1] * msp_factor_i
            for msp_key in msp_signal_labels_dct
        ]
    else:
        msp_signal_mz_lst = []
        msp_signal_i_lst = []
    msp_absent_labels_dct = ms2_info.get("msp_labels", {}).get("absent_labels", {})
    if msp_absent_labels_dct:
        msp_absent_mz_lst = [
            msp_absent_labels_dct[msp_key][0] for msp_key in msp_absent_labels_dct
        ]
        msp_absent_i_lst = [
            msp_absent_labels_dct[msp_key][1] * -1
            # msp_absent_labels_dct[msp_key][1] * msp_factor_i
            for msp_key in msp_absent_labels_dct
        ]
    else:
        msp_absent_mz_lst = []
        msp_absent_i_lst = []
    msp_i_lst = msp_signal_i_lst.copy()
    msp_i_lst.extend(msp_absent_i_lst)

    if msp_i_lst:
        fp_max_i = min(msp_i_lst)
    else:
        fp_max_i = base_peak_i * -0.5

    # plot fingerprint
    min_fp_i = fp_max_i * 1.55  # set min fp y val
    min_fp_h_i = abs(fp_max_i * 0.5)  # set min fp y height
    fingerprint_signal_labels_dct = ms2_info.get("fingerprint_labels", {}).get(
        "signal_labels", {}
    )
    fingerprint_absent_labels_dct = ms2_info.get("fingerprint_labels", {}).get(
        "absent_labels", {}
    )
    if fingerprint_signal_labels_dct:
        fingerprint_signal_mz_lst = [
            fingerprint_signal_labels_dct[fp_key][0]
            for fp_key in fingerprint_signal_labels_dct
        ]
        marker_l, stem_l, base_l = ms2_fig.stem(
            fingerprint_signal_mz_lst,
            [fp_max_i * 1.1] * len(fingerprint_signal_mz_lst),
            linefmt=":",
            markerfmt=" ",
            use_line_collection=True,
        )
        plt.setp(stem_l, color="#00D8D8", alpha=0.4)
        plt.setp(base_l, visible=False)
        for fingerprint_signal_mz in fingerprint_signal_mz_lst:
            obs_fp_box = patches.Rectangle(
                (fingerprint_signal_mz - 1.75, min_fp_i),
                3.5,
                min_fp_h_i,
                facecolor="#00D8D8",
                edgecolor="none",
            )
            ms2_fig.add_patch(obs_fp_box)

    if fingerprint_absent_labels_dct:
        fingerprint_absent_mz_lst = [
            fingerprint_absent_labels_dct[fp_key][0]
            for fp_key in fingerprint_absent_labels_dct
        ]
        marker_l, stem_l, base_l = ms2_fig.stem(
            fingerprint_absent_mz_lst,
            [fp_max_i * 1.1] * len(fingerprint_absent_mz_lst),
            linefmt=":",
            markerfmt=" ",
            use_line_collection=True,
        )
        plt.setp(stem_l, color="#999999", alpha=0.4)
        plt.setp(base_l, visible=False)
        for fingerprint_absent_mz in fingerprint_absent_mz_lst:
            missed_fp_box = patches.Rectangle(
                (fingerprint_absent_mz - 1.75, min_fp_i),
                3.5,
                min_fp_h_i,
                facecolor="#999999",
                edgecolor="none",
            )
            ms2_fig.add_patch(missed_fp_box)

    if msp_signal_mz_lst:
        marker_l, stem_l, base_l = ms2_fig.stem(
            msp_signal_mz_lst,
            msp_signal_i_lst,
            markerfmt=" ",
        )
        plt.setp(stem_l, color="magenta", lw=1, alpha=0.6)
        plt.setp(base_l, visible=False)
    if msp_absent_mz_lst:
        marker_l, stem_l, base_l = ms2_fig.stem(
            msp_absent_mz_lst,
            msp_absent_i_lst,
            markerfmt=" ",
        )
        plt.setp(stem_l, color="grey", lw=1, alpha=0.6)
        plt.setp(base_l, visible=False)

    # plot labels
    txt_props = {"ha": "left", "va": "bottom"}
    observed_specific_peaks = ms2_info.get("scores", {}).get(
        "observed_specific_peaks", {}
    )

    snr_signal_labels_dct = ms2_info.get("snr_labels", {}).get("signal_labels", {})
    # rank_signal_labels_dct = ms2_info.get("rank_labels", {}).get("signal_labels", {})
    # plot_snr_signal_labels_dct.update(rank_signal_labels_dct)
    plot_snr_signal_labels_dct = {}
    for tmp_s_lb in snr_signal_labels_dct:
        if tmp_s_lb.startswith("fp@"):
            pass
        else:
            plot_snr_signal_labels_dct[tmp_s_lb] = snr_signal_labels_dct[tmp_s_lb]
    if plot_snr_signal_labels_dct:
        snr_signal_mz_lst = [
            plot_snr_signal_labels_dct[snr_signal_key][0]
            for snr_signal_key in plot_snr_signal_labels_dct
        ]
        snr_signal_i_lst = [
            plot_snr_signal_labels_dct[snr_signal_key][1]
            for snr_signal_key in plot_snr_signal_labels_dct
        ]
        marker_l, stem_l, base_l = ms2_fig.stem(
            snr_signal_mz_lst, snr_signal_i_lst, markerfmt=" "
        )
        plt.setp(stem_l, color=(0.0, 0.5, 1.0, 0.4), linewidth=1.5)
        plt.setp(base_l, visible=False)
        shift_direction = 1
        for signal_key in plot_snr_signal_labels_dct:
            lb_i = plot_snr_signal_labels_dct[signal_key][1]

            if lb_i/base_peak_i < 0.1:
                lb_i = lb_i + 0.1 * base_peak_i + random.uniform(0, 0.1) * shift_direction * base_peak_i
                shift_direction = shift_direction * -1

            if signal_key.startswith('fp@'):
                lb_color = "#00D8D8"
            else:
                lb_color = (0.0, 0.5, 1.0, 1.0)
            if not signal_key.startswith('fp@'):
                ms2_fig.text(
                    plot_snr_signal_labels_dct[signal_key][0] * 0.973,
                    lb_i,
                    signal_key,
                    txt_props,
                    fontsize=4.5,
                    color=lb_color,
                    rotation=25,
                )
            if signal_key in observed_specific_peaks:
                ms2_fig.text(
                    plot_snr_signal_labels_dct[signal_key][0] - 5,
                    lb_i + 0.224 * base_peak_i,
                    "★",
                    txt_props,
                    fontsize=4.5,
                    color=(0.0, 0.5, 1.0, 1.0),
                )

    snr_noise_labels_dct = ms2_info.get("snr_labels", {}).get("noise_labels", {})
    # rank_noise_labels_dct = ms2_info.get("rank_labels", {}).get("noise_labels", {})
    # plot_snr_noise_labels_dct.update(rank_noise_labels_dct)
    plot_snr_noise_labels_dct = {}
    for tmp_s_lb in snr_noise_labels_dct:
        if tmp_s_lb.startswith("fp@"):
            pass
        else:
            plot_snr_noise_labels_dct[tmp_s_lb] = snr_noise_labels_dct[tmp_s_lb]
    # for tmp_n_lb in plot_snr_noise_labels_dct:
    #     if tmp_n_lb.startswith("fp@"):
    #         plot_snr_noise_labels_dct.pop(tmp_n_lb, None)
    if plot_snr_noise_labels_dct:
        snr_noise_mz_lst = [
            plot_snr_noise_labels_dct[snr_noise_key][0]
            for snr_noise_key in plot_snr_noise_labels_dct
        ]
        snr_noise_i_lst = [
            plot_snr_noise_labels_dct[snr_noise_key][1]
            for snr_noise_key in plot_snr_noise_labels_dct
        ]
        marker_l, stem_l, base_l = ms2_fig.stem(
            snr_noise_mz_lst, snr_noise_i_lst, markerfmt=" "
        )
        plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.4), linewidth=1.5)
        plt.setp(base_l, visible=False)

    # ms2 spectrum
    if not ms2_df.empty:
        ms2_fig.set_xlim(
            [min(ms2_df["mz"].min(), 100) - 50, max(ms2_df["mz"].max() + 50, 800)]
        )
        ms2_fig.set_ylim(
            [min([ms2_df["i"].max() * -1.5, fp_max_i]), ms2_df["i"].max() * 2.25]
        )
        marker_l, stem_l, base_l = ms2_fig.stem(
            ms2_df["mz"].values.tolist(),
            ms2_df["i"].values.tolist(),
            linefmt="black",
            markerfmt=" ",
            basefmt="k-",
            use_line_collection=True,
        )
        plt.setp(stem_l, color="grey", lw=0.75, alpha=1)
        plt.setp(base_l, visible=False)
        ms2_fig.axhline(y=0, color="k", linestyle="-", lw=0.75)

    ms2_fig.set_title(ms2_title_str, color=(0.0, 0.4, 1.0, 1.0), fontsize=8, y=0.98)


def plot_ms2_zoomed_low_fig(
    ms2_zoomed_low_fig, ms2_info: dict, ms2_zoomed_low_title_str: str
):
    # print('start to plot MS/MS <= m/z 400 ...')
    ms2_zoomed_low_fig.tick_params(axis="both", which="major", labelsize=10)
    ms2_zoomed_low_fig.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    ms2_zoomed_low_fig.set_xlabel("m/z", fontsize=7, labelpad=-1)
    ms2_zoomed_low_fig.set_ylabel("Intensity", fontsize=7)
    ms2_zoomed_low_fig.set_title(
        ms2_zoomed_low_title_str, color=(0.0, 0.4, 1.0, 1.0), fontsize=8, y=0.98
    )

    ms2_df = ms2_info.get("df", pd.DataFrame())
    ms2_zoomed_low_df = ms2_df.query(f"mz <= 400")
    txt_props = {"ha": "left", "va": "bottom"}

    pre_snr_signal_labels_dct = ms2_info.get("snr_labels", {}).get("signal_labels", {})
    snr_signal_labels_dct = {}
    for s_key in pre_snr_signal_labels_dct:
        if pre_snr_signal_labels_dct[s_key][0] <= 400:
            snr_signal_labels_dct[s_key] = pre_snr_signal_labels_dct[s_key]
    if snr_signal_labels_dct:
        snr_signal_mz_lst = [
            snr_signal_labels_dct[snr_signal_key][0]
            for snr_signal_key in snr_signal_labels_dct
        ]
        snr_signal_i_lst = [
            snr_signal_labels_dct[snr_signal_key][1]
            for snr_signal_key in snr_signal_labels_dct
        ]
        marker_l, stem_l, base_l = ms2_zoomed_low_fig.stem(
            snr_signal_mz_lst, snr_signal_i_lst, markerfmt=" "
        )
        plt.setp(stem_l, color=(0.0, 0.5, 1.0, 0.4), linewidth=1.5)
        plt.setp(base_l, visible=False)
        for signal_key in snr_signal_labels_dct:
            if signal_key.startswith('fp@'):
                lb_color = "#00D8D8"
            else:
                lb_color = (0.0, 0.5, 1.0, 1.0)

            if not signal_key.startswith('fp@'):
                ms2_zoomed_low_fig.text(
                    snr_signal_labels_dct[signal_key][0] - 5,
                    snr_signal_labels_dct[signal_key][1],
                    signal_key,
                    txt_props,
                    fontsize=4.5,
                    color=lb_color,
                    rotation=60,
                )
                ms2_zoomed_low_fig.text(
                    snr_signal_labels_dct[signal_key][0],
                    snr_signal_labels_dct[signal_key][1],
                    f"{snr_signal_labels_dct[signal_key][0]:.4f}",
                    txt_props,
                    fontsize=4.5,
                    color=lb_color,
                    rotation=60,
                )

    pre_snr_noise_labels_dct = ms2_info.get("snr_labels", {}).get("noise_labels", {})
    snr_noise_labels_dct = {}
    for n_key in pre_snr_noise_labels_dct:
        if pre_snr_noise_labels_dct[n_key][0] <= 400:
            snr_noise_labels_dct[n_key] = pre_snr_noise_labels_dct[n_key]
    if snr_noise_labels_dct:
        snr_noise_mz_lst = [
            snr_noise_labels_dct[snr_noise_key][0]
            for snr_noise_key in snr_noise_labels_dct
        ]
        snr_noise_i_lst = [
            snr_noise_labels_dct[snr_noise_key][1]
            for snr_noise_key in snr_noise_labels_dct
        ]
        marker_l, stem_l, base_l = ms2_zoomed_low_fig.stem(
            snr_noise_mz_lst, snr_noise_i_lst, markerfmt=" "
        )
        plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.4), linewidth=1.5)
        plt.setp(base_l, visible=False)
        for noise_key in snr_noise_labels_dct:
            ms2_zoomed_low_fig.text(
                snr_noise_labels_dct[noise_key][0],
                snr_noise_labels_dct[noise_key][1] * 1.12,
                f"{snr_noise_labels_dct[noise_key][0]:.4f} {noise_key}",
                txt_props,
                fontsize=4.5,
                color="red",
                rotation=90,
            )
            # if "#" in noise_key:
            #     ms2_zoomed_low_fig.text(
            #         snr_noise_labels_dct[noise_key][0] - 5,
            #         snr_noise_labels_dct[noise_key][1],
            #         noise_key,
            #         txt_props,
            #         fontsize=4.5,
            #         color="red",
            #         rotation=90,
            #     )

    # msms spectrum zoomed <= 400 start
    if not ms2_zoomed_low_df.empty:
        ms2_zoomed_low_fig.set_xlim([min(ms2_zoomed_low_df["mz"].min(), 100) - 1, 400])
        ms2_zoomed_low_fig.set_ylim([0, ms2_zoomed_low_df["i"].max() * 2])
        marker_l, stem_l, base_l = ms2_zoomed_low_fig.stem(
            ms2_zoomed_low_df["mz"].values.tolist(),
            ms2_zoomed_low_df["i"].values.tolist(),
            markerfmt=" ",
        )
        plt.setp(stem_l, color="grey", lw=0.75, alpha=0.6)
        plt.setp(base_l, visible=False)


def plot_ms2_zoomed_high_fig(
    ms2_zoomed_high_fig, ms2_info: dict, ms2_zoomed_high_title_str: str
):

    # print('start to plot MS/MS > m/z 400 ...')
    ms2_zoomed_high_fig.tick_params(axis="both", which="major", labelsize=10)
    ms2_zoomed_high_fig.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    ms2_zoomed_high_fig.set_xlabel("m/z", fontsize=7, labelpad=-1)
    ms2_zoomed_high_fig.set_ylabel("Intensity", fontsize=7)
    ms2_zoomed_high_fig.set_title(
        ms2_zoomed_high_title_str, color=(0.0, 0.4, 1.0, 1.0), fontsize=8, y=0.98
    )

    ms2_df = ms2_info.get("df", pd.DataFrame())
    ms2_zoomed_high_df = ms2_df.query(f"mz > 400")

    txt_props = {"ha": "left", "va": "bottom"}

    pre_snr_signal_labels_dct = ms2_info.get("snr_labels", {}).get("signal_labels", {})
    snr_signal_labels_dct = {}
    for s_key in pre_snr_signal_labels_dct:
        if pre_snr_signal_labels_dct[s_key][0] > 400:
            snr_signal_labels_dct[s_key] = pre_snr_signal_labels_dct[s_key]
    if snr_signal_labels_dct:
        snr_signal_mz_lst = [
            snr_signal_labels_dct[snr_signal_key][0]
            for snr_signal_key in snr_signal_labels_dct
        ]
        snr_signal_i_lst = [
            snr_signal_labels_dct[snr_signal_key][1]
            for snr_signal_key in snr_signal_labels_dct
        ]
        marker_l, stem_l, base_l = ms2_zoomed_high_fig.stem(
            snr_signal_mz_lst, snr_signal_i_lst, markerfmt=" "
        )
        plt.setp(stem_l, color=(0.0, 0.5, 1.0, 0.4), linewidth=1.5)
        plt.setp(base_l, visible=False)
        for signal_key in snr_signal_labels_dct:
            if signal_key.startswith('fp@'):
                lb_color = "#00D8D8"
            else:
                lb_color = (0.0, 0.5, 1.0, 1.0)

            if not signal_key.startswith('fp@'):
                ms2_zoomed_high_fig.text(
                    snr_signal_labels_dct[signal_key][0] - 7.5,
                    snr_signal_labels_dct[signal_key][1],
                    signal_key,
                    txt_props,
                    fontsize=4.5,
                    color=lb_color,
                    rotation=60,
                )
                ms2_zoomed_high_fig.text(
                    snr_signal_labels_dct[signal_key][0],
                    snr_signal_labels_dct[signal_key][1],
                    f"{snr_signal_labels_dct[signal_key][0]:.4f}",
                    txt_props,
                    fontsize=4.5,
                    color=lb_color,
                    rotation=60,
                )

    pre_snr_noise_labels_dct = ms2_info.get("snr_labels", {}).get("noise_labels", {})
    snr_noise_labels_dct = {}
    for n_key in pre_snr_noise_labels_dct:
        if pre_snr_noise_labels_dct[n_key][0] > 400:
            snr_noise_labels_dct[n_key] = pre_snr_noise_labels_dct[n_key]
    if snr_noise_labels_dct:
        snr_noise_mz_lst = [
            snr_noise_labels_dct[snr_noise_key][0]
            for snr_noise_key in snr_noise_labels_dct
        ]
        snr_noise_i_lst = [
            snr_noise_labels_dct[snr_noise_key][1]
            for snr_noise_key in snr_noise_labels_dct
        ]
        marker_l, stem_l, base_l = ms2_zoomed_high_fig.stem(
            snr_noise_mz_lst, snr_noise_i_lst, markerfmt=" "
        )
        plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.4), linewidth=1.5)
        plt.setp(base_l, visible=False)
        for noise_key in snr_noise_labels_dct:
            ms2_zoomed_high_fig.text(
                snr_noise_labels_dct[noise_key][0],
                snr_noise_labels_dct[noise_key][1] * 1.12,
                f"{snr_noise_labels_dct[noise_key][0]:.4f} {noise_key}",
                txt_props,
                fontsize=4.5,
                color="red",
                rotation=90,
            )
            # if "#" in noise_key:
            #     ms2_zoomed_high_fig.text(
            #         snr_noise_labels_dct[noise_key][0] - 7.5,
            #         snr_noise_labels_dct[noise_key][1],
            #         noise_key,
            #         txt_props,
            #         fontsize=4.5,
            #         color="red",
            #         rotation=90,
            #     )

    # msms spectrum zoomed > 400 start
    if not ms2_zoomed_high_df.empty:
        ms2_zoomed_high_fig.set_xlim([400, max(ms2_df["mz"].max() + 5, 800)])
        ms2_zoomed_high_fig.set_ylim([0, ms2_zoomed_high_df["i"].max() * 2])
        marker_l, stem_l, base_l = ms2_zoomed_high_fig.stem(
            ms2_zoomed_high_df["mz"].values.tolist(),
            ms2_zoomed_high_df["i"].values.tolist(),
            markerfmt=" ",
        )
        plt.setp(stem_l, color="grey", lw=0.75, alpha=0.6)
        plt.setp(base_l, visible=False)


def plot_spectra(
    figure_info: dict,
    output_path: str,
    dpi=300,
):
    figure_name = figure_info.get("figure", "")
    figure_path = os.path.join(os.path.abspath(output_path), figure_name)

    error_lst = []

    # get lipid info
    info_lipid_dct = figure_info.get("lipid", {})
    abbr = info_lipid_dct.get("lipid_name", "")
    formula_charged = info_lipid_dct.get("charged_formula", "")
    adduct = info_lipid_dct.get("adduct", "")
    lib_mz = info_lipid_dct.get("mz", -1)
    ms1_obs = info_lipid_dct.get("ms1_observed_mz", -1.0)
    ms1_obs_i = info_lipid_dct.get("ms1_observed_i", -1.0)
    ms1_obs_mz = ms1_obs
    ms1_obs_ppm = info_lipid_dct.get("ms1_observed_ppm", -1.0)
    ms1_tolerance_ppm = info_lipid_dct.get("ms1_tolerance_ppm", -1.0)
    ms2_tolerance_ppm = info_lipid_dct.get("ms2_tolerance_ppm", -1.0)

    # get scores info
    info_scores_dct = figure_info.get("scores", {})

    # get spectra info
    info_spectra_dct = figure_info.get("spectra", {})

    xic_dct = info_spectra_dct.get("xic", {})

    ms1_dct = info_spectra_dct.get("ms1", {})
    ms1_scan_id = ms1_dct.get("scan_id", 0)
    ms1_scan_time = ms1_dct.get("scan_time", 0)
    ms1_peaks_dct = ms1_dct.get("peaks", {})
    if ms1_peaks_dct:
        ms1_df = pd.DataFrame(data=ms1_peaks_dct, columns=["mz", "i"])
    else:
        ms1_df = pd.DataFrame()
        error_lst.append("No MS1 spectrum")
    ms1_dct["df"] = ms1_df

    ms2_dct = info_spectra_dct.get("ms2", {})
    ms2_scan_id = ms2_dct.get("scan_id", 0)
    ms2_scan_time = ms2_dct.get("scan_time", 0)
    ms2_pr_mz = ms2_dct.get("precursor_mz", 0)
    ms2_peaks_dct = ms2_dct.get("peaks", {})
    if ms2_peaks_dct:
        ms2_df = pd.DataFrame(data=ms2_peaks_dct, columns=["mz", "i"])
    else:
        ms2_df = pd.DataFrame()
        error_lst.append("No MS2 spectrum")
    ms2_dct["df"] = ms2_df

    # get labels info
    labels_dct = figure_info.get("labels", {})
    isotope_labels_dct = labels_dct.get("isotope_labels", {})
    rank_labels_dct = labels_dct.get("rank_labels", {})
    msp_labels_dct = labels_dct.get("msp_labels", {})
    fingerprint_labels_dct = labels_dct.get("fingerprint_labels", {})
    snr_labels_dct = labels_dct.get("snr_labels", {})

    ms1_dct["isotope_labels"] = isotope_labels_dct
    ms2_dct["rank_labels"] = rank_labels_dct
    ms2_dct["msp_labels"] = msp_labels_dct
    ms2_dct["fingerprint_labels"] = fingerprint_labels_dct
    ms2_dct["snr_labels"] = snr_labels_dct

    # prepare ms1 and ms2 plot data
    ms1_dct["lib_mz"] = lib_mz
    ms1_dct["obs_mz"] = ms1_obs_mz
    ms1_dct["obs_i"] = ms1_obs_i
    ms1_dct["obs_ppm"] = ms1_obs_ppm
    ms1_dct["ms1_tolerance_ppm"] = ms1_tolerance_ppm
    ms1_dct["isotope_score"] = info_scores_dct.get("isotope_score", 0)
    ms1_dct["m'_isotope_score"] = info_scores_dct.get("m'_isotope_score", 0)

    ms2_dct["lipid_name"] = abbr
    ms2_dct["precursor_mz"] = ms2_pr_mz
    ms2_dct["scan_time"] = ms2_scan_time
    ms2_dct["scan_id"] = ms2_scan_time
    ms2_dct["scores"] = info_scores_dct

    if not error_lst:
        # Generate A4 image in landscape
        fig, fig_array = plt.subplots(
            nrows=3, ncols=2, figsize=(11.692, 8.267), sharex="none", sharey="none"
        )

        # init 6 sub figures
        xic_fig = fig_array[0, 0]
        ms1_fig = fig_array[1, 0]
        ms1_zoomed_fig = fig_array[2, 0]
        ms2_fig = fig_array[0, 1]
        ms2_zoomed_low_fig = fig_array[1, 1]
        ms2_zoomed_high_fig = fig_array[2, 1]

        # Make better spacing between subplots
        plt.tight_layout()

        # plot sub figures
        xic_title_str = f"XIC of {abbr} # {adduct} @ lib m/z {lib_mz:.4f}"
        plot_xic(xic_fig, xic_dct, xic_title_str, ms1_scan_time, ms2_scan_time)

        ms1_title_str = f"MS @ {ms1_scan_time:.3f} min | Scan ID: {ms1_scan_id}"
        plot_ms1_fig(ms1_fig, ms1_dct, ms1_title_str, ms1_obs_mz, ms1_obs_i)

        ms_zoomed_title_str = f"Observed m/z {ms1_obs_mz:.4f} ~ {formula_charged} # {adduct} @ {ms1_obs_ppm:.2f}ppm "
        plot_ms1_zoomed_fig(ms1_zoomed_fig, ms1_dct, ms_zoomed_title_str)

        ms2_title_str = f"MS2 of m/z {ms2_pr_mz:.4f} @ {ms2_scan_time:.3f} min | Scan ID: {ms2_scan_id}"
        plot_ms2_fig(ms2_fig, ms2_dct, ms2_title_str)

        plot_ms2_zoomed_low_fig(ms2_zoomed_low_fig, ms2_dct, "MS2 zoomed ≤ m/z 400")
        plot_ms2_zoomed_high_fig(ms2_zoomed_high_fig, ms2_dct, "MS2 zoomed > m/z 400")

        # save to image file
        plt.savefig(figure_path, dpi=dpi)
        print(f"[INFO] --> Image saved as: {figure_path}")
        plt.close()
    else:
        pass

    if os.path.isfile(figure_path):
        has_output_fig = True
    else:
        has_output_fig = False

    return has_output_fig


if __name__ == "__main__":
    import json

    test_folder = r"/Users/lmai/PycharmProjects/tiger1.6/temp"
    js_file_path = r"/Users/lmai/PycharmProjects/tiger1.6/temp/score_info.json"
    with open(js_file_path, "r") as js_obj:
        test_dct = json.load(js_obj)
        plot_spectra(test_dct, test_folder)
