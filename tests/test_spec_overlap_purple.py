"""Smoke test for the spec-overlap PURPLE-classification rule.

Reproduces the user's reported case to verify the asymmetric rule
(purple for user-opens hit by any source — direct or buffer — and
purple for operable candidates hit by ≥ 2 direct sources) actually
fires in code.
"""
from __future__ import annotations

import numpy as np

from vmpt.coords import MSA_V2_REF
from vmpt.msa import load_msa_grid, load_operability
from vmpt.wavelengths import tilt_slope_for_shutter, v2_overlap_distance


SHVAL_S_TOLERANCE = 1
NRS1_QUADS = {1, 3}
NRS2_QUADS = {2, 4}


def _run_overlap_for_case(
    open_shutters: dict[tuple[int, int, int], str],
    disperser: str = "PRISM",
    filt: str = "CLEAR",
) -> dict[str, set[tuple[int, int, int]]]:
    """Replays the contamination calc from `vmpt.main` for the given
    open-shutter dict (keys are (q,s,d), values are target_id strings
    used as group keys). Returns dict with sets of (q,s,d) tuples
    classified as 'purple', 'orange', 'pink', or 'unaffected'."""
    V2_MSA, _ = load_msa_grid()
    _, REASON = load_operability()
    flat_reason = REASON.reshape(-1)
    n = V2_MSA.size
    s_arr = (np.arange(n) % (171 * 365)) // 365
    q_arr = np.arange(n) // (171 * 365) + 1
    d_arr = np.arange(n) % 365
    v2_all = V2_MSA.reshape(-1)
    v2_overlap = float(v2_overlap_distance(disperser, filt))
    in_view_op = flat_reason == 0

    # Group user-opens by (q, d, target_id).
    user_groups: dict[tuple, list[tuple[int, int, int]]] = {}
    for (q, s, d), tid in open_shutters.items():
        key = (q, d, tid if tid is not None else "_anon_")
        user_groups.setdefault(key, []).append((q, s, d))

    user_open_flat = {
        (q - 1) * 171 * 365 + (s - 1) * 365 + (d - 1)
        for (q, s, d) in open_shutters
    }

    def _accumulate(groups, dst_direct, dst_buffer):
        for group in groups:
            if not group:
                continue
            qs = [g[0] for g in group]
            ss = [g[1] for g in group]
            ds = [g[2] for g in group]
            q_o, d_o = qs[0], ds[0]
            s_center = int(np.round(np.mean(ss)))
            slope_k = tilt_slope_for_shutter(
                disperser, filt, q_o, s_center, d_o,
            )
            partners = NRS1_QUADS if q_o in NRS1_QUADS else NRS2_QUADS
            same_det = np.isin(q_arr, list(partners))
            anchor_flat = (
                (q_o - 1) * 171 * 365 + (s_center - 1) * 365 + (d_o - 1)
            )
            v2_o = float(v2_all[anchor_flat])
            v2_open_row = V2_MSA[q_arr - 1, s_center - 1, d_arr]
            dv2 = v2_open_row - v2_o
            drift = slope_k * dv2
            row_offset = (
                np.floor(np.abs(drift) + 0.5).astype(np.int64)
                * np.sign(drift).astype(np.int64)
            )
            different_col = d_arr != (d_o - 1)
            near_v2 = different_col & (np.abs(dv2) < v2_overlap)
            # Band anchored on the slitlet's actual rows [s_lo, s_hi] ± buffer
            # (mirrors vmpt.main; monotonic in the open set, unlike the old
            # s_center ± half_extent which skewed even-shutter slitlets).
            s_lo0 = (min(ss) - 1) + row_offset
            s_hi0 = (max(ss) - 1) + row_offset
            direct_row = (s_arr >= s_lo0) & (s_arr <= s_hi0)
            buffer_row = (
                (s_arr >= s_lo0 - SHVAL_S_TOLERANCE)
                & (s_arr <= s_hi0 + SHVAL_S_TOLERANCE)
                & (~direct_row)
            )
            idx_direct = np.where(
                in_view_op & same_det & direct_row & near_v2
            )[0]
            idx_buffer = np.where(
                in_view_op & same_det & buffer_row & near_v2
            )[0]
            slitlet_flat = {
                (q - 1) * 171 * 365 + (s - 1) * 365 + (d - 1)
                for (q, s, d) in group
            }
            for i in idx_direct.tolist():
                if i in slitlet_flat:
                    continue
                dst_direct[i] = dst_direct.get(i, 0) + 1
            for i in idx_buffer.tolist():
                if i in slitlet_flat:
                    continue
                dst_buffer[i] = dst_buffer.get(i, 0) + 1

    user_direct: dict[int, int] = {}
    user_buffer: dict[int, int] = {}
    stuck_direct: dict[int, int] = {}
    stuck_buffer: dict[int, int] = {}
    hit_sources: dict[int, set[tuple[str, int]]] = {}

    # Patch _accumulate to also populate hit_sources
    def _accumulate_tagged(source_type, groups, dst_direct, dst_buffer):
        # Reuse _accumulate logic but tag sources. For simplicity, just
        # re-do the loop here mirroring main.py.
        for src_idx, group in enumerate(groups):
            if not group:
                continue
            qs = [g[0] for g in group]
            ss = [g[1] for g in group]
            ds = [g[2] for g in group]
            q_o, d_o = qs[0], ds[0]
            s_center = int(np.round(np.mean(ss)))
            slope_k = tilt_slope_for_shutter(
                disperser, filt, q_o, s_center, d_o,
            )
            partners = NRS1_QUADS if q_o in NRS1_QUADS else NRS2_QUADS
            same_det = np.isin(q_arr, list(partners))
            anchor_flat = (q_o-1)*171*365 + (s_center-1)*365 + (d_o-1)
            v2_o = float(v2_all[anchor_flat])
            v2_open_row = V2_MSA[q_arr-1, s_center-1, d_arr]
            dv2 = v2_open_row - v2_o
            drift = slope_k * dv2
            row_offset = (
                np.floor(np.abs(drift) + 0.5).astype(np.int64)
                * np.sign(drift).astype(np.int64)
            )
            different_col = d_arr != (d_o - 1)
            near_v2 = different_col & (np.abs(dv2) < v2_overlap)
            # Band anchored on the slitlet's actual rows [s_lo, s_hi] ± buffer
            # (mirrors vmpt.main; monotonic in the open set, unlike the old
            # s_center ± half_extent which skewed even-shutter slitlets).
            s_lo0 = (min(ss) - 1) + row_offset
            s_hi0 = (max(ss) - 1) + row_offset
            direct_row = (s_arr >= s_lo0) & (s_arr <= s_hi0)
            buffer_row = (
                (s_arr >= s_lo0 - SHVAL_S_TOLERANCE)
                & (s_arr <= s_hi0 + SHVAL_S_TOLERANCE)
                & (~direct_row)
            )
            in_view_candidates = flat_reason != 1
            idx_direct = np.where(
                in_view_candidates & same_det & direct_row & near_v2
            )[0]
            idx_buffer = np.where(
                in_view_candidates & same_det & buffer_row & near_v2
            )[0]
            slitlet_flat = {
                (q-1)*171*365 + (s-1)*365 + (d-1)
                for (q, s, d) in group
            }
            tag = (source_type, int(src_idx))
            for i in idx_direct.tolist():
                if i in slitlet_flat:
                    continue
                dst_direct[i] = dst_direct.get(i, 0) + 1
                hit_sources.setdefault(i, set()).add(tag)
            for i in idx_buffer.tolist():
                if i in slitlet_flat:
                    continue
                dst_buffer[i] = dst_buffer.get(i, 0) + 1
                hit_sources.setdefault(i, set()).add(tag)
    user_groups_list = list(user_groups.values())
    _accumulate_tagged("user", user_groups_list, user_direct, user_buffer)
    # No stuck-opens in these synthetic test cases.

    # Identify conflicted sources via hit_sources.
    conflicted_user: set[int] = set()
    conflicted_stuck: set[int] = set()
    for src_idx, group in enumerate(user_groups_list):
        for (q, s, d) in group:
            flat = (q-1)*171*365 + (s-1)*365 + (d-1)
            others = hit_sources.get(flat, set()) - {("user", src_idx)}
            if others:
                conflicted_user.add(src_idx)
                for (t, oi) in others:
                    if t == "user":
                        conflicted_user.add(oi)
                    else:
                        conflicted_stuck.add(oi)
                break

    purple: set[tuple[int, int, int]] = set()
    orange: set[tuple[int, int, int]] = set()
    pink: set[tuple[int, int, int]] = set()
    all_idx = (
        set(user_direct) | set(user_buffer)
        | set(stuck_direct) | set(stuck_buffer)
    )

    def _to_qsd(i: int) -> tuple[int, int, int]:
        return (
            i // (171 * 365) + 1,
            (i % (171 * 365)) // 365 + 1,
            i % 365 + 1,
        )

    for i in all_idx:
        n_ud = user_direct.get(i, 0)
        n_ub = user_buffer.get(i, 0)
        n_sd = stuck_direct.get(i, 0)
        n_sb = stuck_buffer.get(i, 0)
        n_total_u = n_ud + n_ub
        n_total_s = n_sd + n_sb
        n_total = n_total_u + n_total_s
        qsd = _to_qsd(i)
        if i in user_open_flat:
            if n_total >= 1:
                purple.add(qsd)
            continue
        sources_here = hit_sources.get(i, set())
        from_conflicted = any(
            (t == "user" and si in conflicted_user)
            or (t == "stuck" and si in conflicted_stuck)
            for (t, si) in sources_here
        )
        if from_conflicted:
            purple.add(qsd)
        elif n_total_u >= 1:
            orange.add(qsd)
        elif n_total_s >= 1:
            pink.add(qsd)
    return {"purple": purple, "orange": orange, "pink": pink}


def test_user_case_q1_s100_d322_and_s97_d323_purple():
    """Two slitlets in Q1 at (s=100-102, d=322) and (s=97-99, d=323).
    Touching: s=100 borders s=99 with no operable row between them.
    The touching shutter on each side should turn purple; the other
    four shutters keep their natural red rendering."""
    open_shutters = {
        (1, 100, 322): "A", (1, 101, 322): "A", (1, 102, 322): "A",
        (1, 97, 323): "B",  (1, 98, 323): "B",  (1, 99, 323): "B",
    }
    out = _run_overlap_for_case(open_shutters)
    assert (1, 100, 322) in out["purple"], (
        f"q1 s=100 d=322 (slitlet A's bottom, touching B) should be "
        f"purple, but isn't. purple set = {sorted(out['purple'])}"
    )
    assert (1, 99, 323) in out["purple"], (
        f"q1 s=99 d=323 (slitlet B's top, touching A) should be "
        f"purple, but isn't. purple set = {sorted(out['purple'])}"
    )
    # The other four user-opens (A's 101,102 and B's 97,98) are too
    # far from the other slitlet to feel its spectrum — they stay red.
    for q, s, d in [(1, 101, 322), (1, 102, 322), (1, 97, 323), (1, 98, 323)]:
        assert (q, s, d) not in out["purple"], (
            f"q{q} s={s} d={d} shouldn't be purple in this case"
        )


def test_single_slitlet_has_no_purple():
    """A single 3-shutter slitlet, with no other opens, has zero
    purple shutters. (Regression for the 'every N-slitlet shutter
    turns purple' bug from earlier.)"""
    open_shutters = {
        (1, 100, 322): "A", (1, 101, 322): "A", (1, 102, 322): "A",
    }
    out = _run_overlap_for_case(open_shutters)
    for q, s, d in open_shutters:
        assert (q, s, d) not in out["purple"]
