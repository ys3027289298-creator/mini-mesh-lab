#!/usr/bin/env python3
"""分块流式加载帧级仿真（带宽共享 + 单解码槽 模型）。

两条判定都写成"剩余时间"的函数：
  预取:  T_enter(c) <= k * v_est * (t_net + t_dec)
  降级:  T_dwell(c) < T_hold
衔接: 预取应早于降级  Δ = t_net + t_dec + T_hold - L/v（见 README）。

与真实流式系统一致的两个瓶颈：
  1) 网络带宽 B 在所有在途(fetching)块间均分；越早发起的块在速度突变前
     已累积传输更多字节，因此预取策略能在加速期拉开差距。
  2) 单解码槽串行处理传输完成的块，按紧急度(T_enter)调度；槽满则排队，
     排队期间块进入视野即记缺块。
块状态: unseen -> fetching -> queued -> loading -> active -> degraded。
"""
import argparse
import json
import math
import os
import random

S = 10.0           # 块边长 (m)
D = 40.0           # 视野深度 (m)
CHUNK_BYTES = 1_000_000      # 每块字节数
BANDWIDTH = 18_000_000       # 总网络带宽 (B/s)，在途块均分
T_DEC = 0.20       # 解码时间基准 (s)
T_DEC_JITTER = 0.06
T_HOLD = 0.60      # 降级安全时间 (s)
FPS = 60
DURATION = 120.0
MAX_INFLIGHT = 1   # 解码槽数
PREFETCH_MARGIN = 1.15
TRIGGER_SLACK = 1.10
SCAN_AHEAD = 28
V_WINDOW = 0.3     # 瞬时策略滚动峰值窗口 (s)

SPEED_TIERS = [
    # 名称, 巡航平均速度(m/s), 加速段目标速度, 振荡幅值
    ("低速", 4.0, 8.0, 0.8),
    ("中速", 12.0, 20.0, 2.0),
    ("高速", 24.0, 38.0, 4.0),
]


def velocity(t, vbar, v_target, amp):
    """20s 周期：巡航 8s -> 2s 急加速 -> 高速巡航 8s -> 2s 减速。"""
    phase = t % 20.0
    if phase < 8.0:
        base = vbar
    elif phase < 10.0:
        base = vbar + (v_target - vbar) * ((phase - 8.0) / 2.0)
    elif phase < 18.0:
        base = v_target
    else:
        base = v_target - (v_target - vbar) * ((phase - 18.0) / 2.0)
    return max(0.5, base + amp * 0.25 * math.sin(2 * math.pi * t / 1.3))


def run_tier(vbar, v_target, amp, policy, seed):
    rng = random.Random(seed)
    dt = 1.0 / FPS
    nframes = int(DURATION * FPS)
    x = 0.0

    state = {}
    fetched_bytes = {}   # 块已传输字节
    ready_at = {}
    prefetch_t = {}
    lead_at_prefetch = {}
    enter_t = {}
    inflight = {}
    v_history = []

    immediate_degrade = 0
    missing_frames = 0
    v_peak = 0.0
    initial_chunks = set()

    # 冷启动：起步块传输已完成，直接占解码槽
    for cid in range(MAX_INFLIGHT):
        state[cid] = 'loading'
        prefetch_t[cid] = 0.0
        initial_chunks.add(cid)
        inflight[cid] = T_DEC + rng.random() * T_DEC_JITTER

    def t_enter_of(cid):
        return max(0.0, (cid * S - x) / max(v, 1e-6))

    for f in range(nframes):
        t = f * dt
        v = velocity(t, vbar, v_target, amp)
        v_peak = max(v_peak, v)
        v_history = [(tt, vv) for tt, vv in v_history if t - tt <= V_WINDOW]
        v_history.append((t, v))
        x += v * dt
        front = math.floor(x / S)

        # 1) 解码完成
        for cid in list(inflight):
            if inflight[cid] <= t:
                state[cid] = 'active'
                ready_at[cid] = inflight.pop(cid)

        # 2) 带宽均分给在途块，累积传输；传完入解码等待队列
        fetching = [c for c, st in state.items() if st == 'fetching']
        if fetching:
            share = BANDWIDTH / len(fetching) * dt
            for cid in fetching:
                fetched_bytes[cid] += share
                if fetched_bytes[cid] >= CHUNK_BYTES:
                    state[cid] = 'queued'

        # 3) 空闲解码槽按紧急度受理 queued 块
        queued = sorted((c for c, st in state.items() if st == 'queued'),
                        key=t_enter_of)
        for cid in queued:
            if len(inflight) >= MAX_INFLIGHT:
                break
            state[cid] = 'loading'
            inflight[cid] = t + T_DEC + rng.random() * T_DEC_JITTER

        # 4) 发起新预取；单块基线传输时间 = CHUNK_BYTES/BANDWIDTH
        t_net = CHUNK_BYTES / BANDWIDTH
        t_load = t_net + T_DEC + 0.5 * T_DEC_JITTER
        if policy == 'avg':
            v_global = x / max(t, dt)
            threshold = v_global * t_load
        else:
            recent_peak = max(vv for _, vv in v_history)
            threshold = TRIGGER_SLACK * PREFETCH_MARGIN * recent_peak * t_load
        scan_horizon = math.floor((x + SCAN_AHEAD * S) / S)
        for cid in range(front + 1, scan_horizon + 1):
            if cid in state:
                continue
            te = t_enter_of(cid)
            if te > threshold:
                continue
            state[cid] = 'fetching'
            fetched_bytes[cid] = 0.0
            prefetch_t[cid] = t
            lead_at_prefetch[cid] = te

        # 5) 缺块：相机已越过块前边缘，块未 active/degraded（冷启动除外）
        for cid in range(max(0, front - int(D / S)), front + 1):
            if cid not in enter_t:
                enter_t[cid] = t
            if cid in initial_chunks:
                continue
            if state.get(cid) not in ('active', 'degraded'):
                missing_frames += 1

        # 6) 降级：预计剩余驻留时间低于 T_hold
        for cid in list(state):
            if state[cid] != 'active':
                continue
            t_dwell = ((cid + 1) * S - x) / max(v, 1e-6)
            if t_dwell < T_HOLD:
                state[cid] = 'degraded'
                if t - ready_at[cid] < T_HOLD:
                    immediate_degrade += 1

    leads = [lead_at_prefetch[c] for c in lead_at_prefetch
             if c in enter_t and c not in initial_chunks]
    lead_mean = sum(leads) / len(leads) if leads else 0.0
    trans = [lead_at_prefetch[c] for c in lead_at_prefetch
             if c in enter_t and c not in initial_chunks
             and 8.0 <= (prefetch_t[c] % 20.0) < 10.0]
    lead_trans = sum(trans) / len(trans) if trans else 0.0
    true_leads = [enter_t[c] - prefetch_t[c] for c in prefetch_t
                  if c in enter_t and c not in initial_chunks]
    true_lead_mean = sum(true_leads) / len(true_leads) if true_leads else 0.0

    return {
        "tier": next(n for n, vb, _, _ in SPEED_TIERS if vb == vbar),
        "policy": policy,
        "vbar": vbar,
        "v_peak": round(v_peak, 2),
        "prefetch_lead_s": round(lead_mean, 3),
        "transition_lead_s": round(lead_trans, 3),
        "true_lead_s": round(true_lead_mean, 3),
        "immediate_degrade_chunks": immediate_degrade,
        "missing_frames": missing_frames,
        "chunks_prefetched": len(prefetch_t),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    rows = []
    for name, vbar, v_target, amp in SPEED_TIERS:
        for policy, seed in (("avg", 101), ("inst", 202)):
            rows.append(run_tier(vbar, v_target, amp, policy, seed))
    text = json.dumps(rows, ensure_ascii=False, indent=2)
    print(text)
    if args.out:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as fh:
            fh.write(text + "\n")


if __name__ == "__main__":
    main()
