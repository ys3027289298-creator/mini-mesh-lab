#!/usr/bin/env python3
"""Render a self-contained interactive HTML report (offline, no dependencies
beyond the stdlib) from the benchmark CSVs and analysis.json.
Usage: python3 -P make_report.py <out_dir>"""
import csv
import json
import sys
from pathlib import Path

OUT = Path(sys.argv[1])

def read_med(pattern, col):
    rows = {}
    for fp in sorted(OUT.glob(pattern)):
        if "trivial" in fp.name:
            continue
        with open(fp, newline="") as f:
            for r in csv.DictReader(f):
                rows.setdefault(int(r["n"]), []).append(float(r[col]))
    ns = sorted(rows)
    return ns, [sorted(rows[n])[len(rows[n]) // 2] / 1e6 for n in ns]

ns, total = read_med("bench_mesh*.csv", "total_ns")
_, proxy = read_med("bench_mesh*.csv", "proxy_ns")
_, broad = read_med("bench_mesh*.csv", "broad_ns")
_, narrow = read_med("bench_mesh*.csv", "narrow_ns")
ext_n, ext_t = read_med("bench_ext*.csv", "total_ns")
ana = json.loads((OUT / "analysis.json").read_text())

series = [
    ("total",  "单帧总耗时", "#52514e", ns, total, False),
    ("proxy",  "代理生成 O(n)", "#2a78d6", ns, proxy, False),
    ("broad",  "粗检测 O(n²)", "#eb6834", ns, broad, False),
    ("narrow", "精检测", "#1baf7a", ns, narrow, True),
]

data_js = json.dumps({
    "ns": ns, "series": [{"id": i, "name": n, "color": c,
                          "x": x, "y": y, "dash": d}
                         for i, (_, n, c, x, y, d) in enumerate(series)],
    "ext": {"x": ext_n, "y": ext_t},
    "cross": ana["proxy_broad_crossover_n_interp"],
    "knee": ana["knee_n_piecewise_fit"],
    "knee15": ana["knee_n_local_exp_1p5"],
    "budget": ana["frame_budget_ms"],
}, ensure_ascii=False)

table_rows = "\n".join(
    f"<tr><td>{n}</td><td>{t:.4f}</td><td>{p:.4f}</td><td>{b:.4f}</td>"
    f"<td>{nr:.5f}</td><td>{100*p/t:.0f}%</td><td>{100*b/t:.0f}%</td></tr>"
    for n, t, p, b, nr in zip(ns, total, proxy, broad, narrow))

html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>碰撞检测代理压测报告</title>
<style>
:root {{
  --surface: #fcfcfb; --panel: #ffffff; --ink: #0b0b0b; --ink2: #52514e;
  --grid: #e3e2de; --band: #d9d8d2; --blue:#2a78d6; --orange:#eb6834;
  --aqua:#1baf7a; --red:#e34948;
}}
@media (prefers-color-scheme: dark) {{
  :root {{ --surface:#1a1a19; --panel:#222221; --ink:#f0efec; --ink2:#c3c2b7;
           --grid:#34342f; --band:#3a3a35; --blue:#3987e5; --orange:#d95926;
           --aqua:#199e70; --red:#e66767; }}
}}
* {{ box-sizing: border-box; }}
body {{ margin:0; background:var(--surface); color:var(--ink);
       font-family:"Noto Sans CJK SC","PingFang SC","Microsoft YaHei",sans-serif; }}
.wrap {{ max-width:1080px; margin:0 auto; padding:28px 24px 60px; }}
h1 {{ font-size:21px; margin:0 0 4px; }}
.sub {{ color:var(--ink2); font-size:13px; margin-bottom:22px; }}
.cards {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(200px,1fr));
          gap:12px; margin-bottom:22px; }}
.card {{ background:var(--panel); border:1px solid var(--grid); border-radius:10px;
         padding:14px 16px; }}
.card .k {{ font-size:12px; color:var(--ink2); margin-bottom:6px; }}
.card .v {{ font-size:22px; font-weight:700; }}
.card .d {{ font-size:12px; color:var(--ink2); margin-top:4px; line-height:1.5; }}
.panel {{ background:var(--panel); border:1px solid var(--grid); border-radius:10px;
          padding:16px 18px; margin-bottom:20px; }}
.panel h2 {{ font-size:15px; margin:0 0 10px; }}
svg {{ width:100%; height:auto; display:block; }}
.legend {{ display:flex; flex-wrap:wrap; gap:16px; font-size:13px; margin-top:8px; }}
.legend span {{ display:inline-flex; align-items:center; gap:6px; color:var(--ink2); }}
.sw {{ width:18px; height:3px; border-radius:2px; }}
#tip {{ position:fixed; pointer-events:none; background:var(--panel);
        border:1px solid var(--grid); border-radius:8px; padding:8px 10px;
        font-size:12px; display:none; z-index:10; box-shadow:0 4px 14px rgba(0,0,0,.18);
        white-space:nowrap; }}
#tip b {{ font-size:13px; }} #tip .r {{ display:flex; gap:8px; align-items:center; }}
#tip i {{ width:10px; height:10px; border-radius:50%; display:inline-block; }}
table {{ width:100%; border-collapse:collapse; font-size:12.5px; }}
th,td {{ padding:6px 8px; text-align:right; border-bottom:1px solid var(--grid); }}
th:first-child,td:first-child {{ text-align:left; }}
th {{ color:var(--ink2); font-weight:600; }}
details {{ margin-top:6px; }} summary {{ cursor:pointer; color:var(--ink2); font-size:13px; }}
.note {{ font-size:13px; line-height:1.75; color:var(--ink2); }}
.note b {{ color:var(--ink); }}
code {{ background:var(--grid); padding:1px 5px; border-radius:4px; font-size:12px; }}
</style>
</head>
<body>
<div class="wrap">
  <h1>碰撞检测代理生成 · 压力测试报告</h1>
  <div class="sub">物体数 10 → 2000（另测 2500–4000 外推点）· 场景密度固定 L=∛n ·
    每点 9 组 × 自适应帧数取中位，3 次独立跑再取中位 · AMD Ryzen 9700X 单核 · gcc -O2</div>

  <div class="cards">
    <div class="card"><div class="k">瓶颈易主拐点</div>
      <div class="v">n ≈ {ana['proxy_broad_crossover_n_interp']}</div>
      <div class="d">粗检测耗时追平代理生成；总阶数在 n={ana['knee_n_local_exp_1p5']} 越过 1.5</div></div>
    <div class="card"><div class="k">拐点前阶数 / 瓶颈</div>
      <div class="v">O(n<sup>1.02</sup>)</div>
      <div class="d">代理生成：256 顶点/物体的旋转拟合 AABB，占帧时间 85–95%</div></div>
    <div class="card"><div class="k">拐点后阶数 / 瓶颈</div>
      <div class="v">O(n<sup>≈2</sup>)</div>
      <div class="d">粗检测 n(n−1)/2 次 AABB 两两测试，n=2000 时占 55%</div></div>
    <div class="card"><div class="k">60 FPS 帧预算上限</div>
      <div class="v">≈ {ana['max_objects_measured_interp']} 个</div>
      <div class="d">实测插值 {ana['max_objects_measured_interp']}，二次拟合 {ana['max_objects_quadratic_fit']}；3500 个 14.2ms，4000 个 17.9ms</div></div>
  </div>

  <div class="panel">
    <h2>单帧耗时 vs 物体数（双对数，悬浮查看数值）</h2>
    <svg id="chart" viewBox="0 0 980 520" role="img" aria-label="log-log timing curve"></svg>
    <div class="legend" id="legend"></div>
  </div>

  <div class="panel">
    <h2>结论</h2>
    <div class="note">
      <p><b>拐点位置：n ≈ 1300–1700 的过渡带，瓶颈易主点 n≈1600。</b>
      分段拟合在 n≈{ana['knee_n_piecewise_fit']} 处首次检测到二次项（该点恰好在本机 L2 缓存
      容量边界：顶点池 3KB/物体 × {ana['knee_n_piecewise_fit']} ≈ 1.35MB，略超 1MB L2）；
      到 n≈{ana['proxy_broad_crossover_n_interp']} 粗检测追平代理生成，总曲线局部阶数在
      n={ana['knee_n_local_exp_1p5']} 越过 1.5。</p>
      <p><b>拐点前（n ≲ 460，阶数 ≈1.0）：</b>瓶颈是<b>代理生成</b>——每帧把每个物体的
      256 个局部顶点做四元数旋转再拟合 AABB，O(n) 但常数大，占帧时间 85–95%，且顶点池
      依次溢出 L1d（n≈16）和 L2（n≈341）。</p>
      <p><b>拐点后（n ≳ 1700，阶数 →≈2）：</b>瓶颈换成<b>粗检测的暴力两两 AABB 测试</b>，
      测试次数 n(n−1)/2 呈二次膨胀（n=2000 时 200 万次/帧，占 55%）；大 n 窗口实测阶数
      收敛到 ~2.1。精检测始终可忽略（稀疏场景下候选对 ≈ 0–5 个/帧，&lt;0.001ms）。</p>
      <p><b>一句话结论：</b>按 16.7ms（60 FPS）帧预算，这个场景<b>最多支持约
      {ana['max_objects_measured_interp']} 个物体</b>（实测插值；保守取整 ≈3800）；
      若要继续扩容，应先把粗检测换成空间划分（网格/排序扫描/SAP），而不是优化代理生成。</p>
    </div>
  </div>

  <div class="panel">
    <details><summary>数据表（中位单帧耗时 ms / 阶段占比）</summary>
    <table>
      <thead><tr><th>n</th><th>总耗时</th><th>代理生成</th><th>粗检测</th>
      <th>精检测</th><th>代理占比</th><th>粗检占比</th></tr></thead>
      <tbody>{table_rows}</tbody>
    </table></details>
  </div>
</div>
<div id="tip"></div>

<script>
const D = {data_js};
const W=980,H=520,M={{l:78,r:30,t:24,b:58}};
const PW=W-M.l-M.r, PH=H-M.t-M.b;
const allx = D.ns.concat(D.ext.x), ally = [];
D.series.forEach(s => ally.push(...s.y)); ally.push(...D.ext.y, D.budget);
const lg10 = x => Math.log10(x);
const xmin=Math.pow(10,Math.floor(lg10(Math.min(...allx)))),
      xmax=Math.pow(10,Math.ceil(lg10(Math.max(...allx)))),
      ymin=Math.pow(10,Math.floor(lg10(Math.min(...ally.filter(v=>v>0))))),
      ymax=Math.pow(10,Math.ceil(lg10(Math.max(...ally))));
const X = v => M.l + (lg10(v)-lg10(xmin))/(lg10(xmax)-lg10(xmin))*PW;
const Y = v => M.t + PH - (lg10(v)-lg10(ymin))/(lg10(ymax)-lg10(ymin))*PH;
const svg = document.getElementById('chart');
const NS='http://www.w3.org/2000/svg';
function el(t,a){{const e=document.createElementNS(NS,t);for(const k in a)e.setAttribute(k,a[k]);return e;}}
// grid + axes
for(let p=ymin;p<=ymax*1.001;p*=10){{
  svg.appendChild(el('line',{{x1:M.l,x2:W-M.r,y1:Y(p),y2:Y(p),stroke:'var(--grid)','stroke-width':1}}));
  const t=el('text',{{x:M.l-8,y:Y(p)+4,'text-anchor':'end','font-size':11,fill:'var(--ink2)'}});
  t.textContent=p>=1?p.toFixed(0):p.toExponential(0); svg.appendChild(t);
}}
for(let p=xmin;p<=xmax*1.001;p*=10){{
  svg.appendChild(el('line',{{y1:M.t,y2:H-M.b,x1:X(p),x2:X(p),stroke:'var(--grid)','stroke-width':1}}));
  const t=el('text',{{x:X(p),y:H-M.b+20,'text-anchor':'middle','font-size':11,fill:'var(--ink2)'}});
  t.textContent=p; svg.appendChild(t);
}}
// transition band
const band=el('rect',{{x:X(D.knee),y:M.t,width:X(D.cross)-X(D.knee),height:PH,fill:'var(--band)',opacity:.45}});
svg.appendChild(band);
// budget line
svg.appendChild(el('line',{{x1:M.l,x2:W-M.r,y1:Y(D.budget),y2:Y(D.budget),
  stroke:'var(--red)','stroke-width':1.6,'stroke-dasharray':'7 5'}}));
const bt=el('text',{{x:M.l+6,y:Y(D.budget)-7,'font-size':12,fill:'var(--red)'}});
bt.textContent='16.7ms 帧预算（60 FPS）'; svg.appendChild(bt);
// handover line
svg.appendChild(el('line',{{x1:X(D.cross),x2:X(D.cross),y1:M.t,y2:H-M.b,
  stroke:'var(--aqua)','stroke-width':1.4,'stroke-dasharray':'4 4'}}));
// series paths + points
D.series.forEach(s=>{{
  const pts=s.x.map((v,i)=>X(v)+','+Y(s.y[i])).join(' ');
  svg.appendChild(el('polyline',{{points:pts,fill:'none',stroke:s.color,'stroke-width':2.2,
    'stroke-dasharray':s.dash?'5 5':'0'}}));
  s.x.forEach((v,i)=>svg.appendChild(el('circle',{{cx:X(v),cy:Y(s.y[i]),r:4.5,fill:s.color,
    stroke:'var(--panel)','stroke-width':1.5}})));
}});
// extrapolation points
D.ext.x.forEach((v,i)=>svg.appendChild(el('circle',{{cx:X(v),cy:Y(D.ext.y[i]),r:6,
  fill:'var(--panel)',stroke:'var(--ink)','stroke-width':2}})));
// axis titles
const tx=el('text',{{x:M.l+PW/2,y:H-12,'text-anchor':'middle','font-size':13,fill:'var(--ink2)'}});
tx.textContent='物体数 n（对数轴）'; svg.appendChild(tx);
const ty=el('text',{{x:18,y:M.t+PH/2,'text-anchor':'middle','font-size':13,fill:'var(--ink2)',
  transform:'rotate(-90 18 '+ (M.t+PH/2) +')'}});
ty.textContent='单帧耗时 ms（对数轴）'; svg.appendChild(ty);
// legend
const lg=document.getElementById('legend');
D.series.forEach(s=>{{
  const sp=document.createElement('span');
  sp.innerHTML='<i class="sw" style="background:'+s.color+';'+(s.dash?'border-top:2px dashed':'')+'"></i>'+s.name;
  lg.appendChild(sp);
}});
lg.insertAdjacentHTML('beforeend','<span><i class="sw" style="background:var(--ink);border-radius:50%;width:10px;height:10px"></i>外推实测点（&gt;2000）</span>');
// crosshair + tooltip
const tip=document.getElementById('tip');
const cross=el('line',{{y1:M.t,y2:H-M.b,stroke:'var(--ink2)','stroke-width':1,
  'stroke-dasharray':'3 3',opacity:0}}); svg.appendChild(cross);
svg.addEventListener('mousemove',e=>{{
  const r=svg.getBoundingClientRect();
  const mx=(e.clientX-r.left)/r.width*W;
  let best=0,bd=1e9;
  D.ns.forEach((v,i)=>{{const d=Math.abs(X(v)-mx);if(d<bd){{bd=d;best=i;}}}});
  const n=D.ns[best];
  cross.setAttribute('x1',X(n));cross.setAttribute('x2',X(n));cross.setAttribute('opacity',.7);
  let html='<b>n = '+n+'</b><br>';
  D.series.forEach(s=>{{html+='<div class="r"><i style="background:'+s.color+'"></i>'+
    s.name+'：'+s.y[best].toFixed(5)+' ms</div>';}});
  tip.innerHTML=html; tip.style.display='block';
  tip.style.left=Math.min(e.clientX+14,innerWidth-tip.offsetWidth-12)+'px';
  tip.style.top=(e.clientY+14)+'px';
}});
svg.addEventListener('mouseleave',()=>{{tip.style.display='none';cross.setAttribute('opacity',0);}});
</script>
</body>
</html>"""

(OUT / "report.html").write_text(html, encoding="utf-8")
print("wrote", OUT / "report.html")
