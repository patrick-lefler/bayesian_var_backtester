(() => {
  'use strict';
  const DAYS = 250, WIN = 25, NW = 10, CLAIM = 0.01;
  const DEFAULTS = {
    b1: 3, d1: 250, b2: 2, d2: 250, w1: 0.4, w2: 0.5, mode: 'history', claimN: 500,
    wins: [0, 1, 0, 0, 0, 3, 3, 0, 0, 0]
  };
  // 1996 Basel traffic light: multiplier by year-end exception count (250 days, 99% VaR)
  const MULT = n => n <= 4 ? 3.00 : n >= 10 ? 4.00 : [3.40, 3.50, 3.65, 3.75, 3.85][n - 5];
  const zoneOf = n => n <= 4 ? 'g' : n <= 9 ? 'y' : 'r';
  const ZNAME = { g: 'Green', y: 'Yellow', r: 'Red' };
  const $ = id => document.getElementById(id);
  const pc = (v, d = 2) => isFinite(v) ? (v * 100).toFixed(d) + '%' : '—';
  const pc0 = v => Math.round(v * 100) + '%';
  const f1 = v => v.toFixed(1), f2 = v => v.toFixed(2);

  // ---------- math ----------
  const LG = [676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61502916214059,12.507343278686905,-0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7];
  function lgamma(x){
    if (x < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * x)) - lgamma(1 - x);
    x -= 1; let a = 0.99999999999980993; const t = x + 7.5;
    for (let i = 0; i < 8; i++) a += LG[i] / (x + i + 1);
    return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
  }
  const lbeta = (a, b) => lgamma(a) + lgamma(b) - lgamma(a + b);
  const lchoose = (n, k) => lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1);
  function betaBinomPmf(R, a, b){
    const out = new Array(R + 1), lb = lbeta(a, b);
    for (let y = 0; y <= R; y++) out[y] = Math.exp(lchoose(R, y) + lbeta(y + a, R - y + b) - lb);
    return out;
  }
  // Regularized incomplete beta I_x(a,b) (continued fraction, Numerical Recipes)
  function betacf(x, a, b){
    const FPMIN = 1e-300; let qab = a + b, qap = a + 1, qam = a - 1, c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < FPMIN) d = FPMIN; d = 1 / d; let h = d;
    for (let m = 1; m <= 300; m++){
      const m2 = 2 * m; let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN; c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN; d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN; c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN; d = 1 / d;
      const del = d * c; h *= del; if (Math.abs(del - 1) < 1e-12) break;
    }
    return h;
  }
  function betaCdf(x, a, b){
    if (x <= 0) return 0; if (x >= 1) return 1;
    const bt = Math.exp(a * Math.log(x) + b * Math.log(1 - x) - lbeta(a, b));
    return x < (a + 1) / (a + b + 2) ? bt * betacf(x, a, b) / a : 1 - bt * betacf(1 - x, b, a) / b;
  }
  function betaQuantile(p, a, b){
    let lo = 0, hi = 1;
    for (let i = 0; i < 80; i++){ const mid = (lo + hi) / 2; if (betaCdf(mid, a, b) < p) lo = mid; else hi = mid; }
    return (lo + hi) / 2;
  }
  // Upper tail of chi-square(k) via regularized gamma Q(k/2, x/2)
  function gammaQ(a, x){
    if (x <= 0) return 1;
    const gln = lgamma(a);
    if (x < a + 1){
      let sum = 1 / a, del = sum, ap = a;
      for (let n = 0; n < 500; n++){ ap += 1; del *= x / ap; sum += del; if (Math.abs(del) < Math.abs(sum) * 1e-14) break; }
      return 1 - sum * Math.exp(-x + a * Math.log(x) - gln);
    }
    let b = x + 1 - a, c = 1e300, d = 1 / b, h = d;
    for (let i = 1; i < 500; i++){
      const an = -i * (i - a); b += 2; d = an * d + b; if (Math.abs(d) < 1e-300) d = 1e-300;
      c = b + an / c; if (Math.abs(c) < 1e-300) c = 1e-300; d = 1 / d; const del = d * c; h *= del; if (Math.abs(del - 1) < 1e-14) break;
    }
    return Math.exp(-x + a * Math.log(x) - gln) * h;
  }
  const chi2Upper = (x, k) => gammaQ(k / 2, x / 2);

  // Dispersion test across 25-day windows observed so far
  function dispersion(wins){
    const k = wins.length, tot = wins.reduce((s, x) => s + x, 0);
    if (k < 3 || tot < 3) return { ok: false, k, tot };
    const ph = tot / (WIN * k), e = WIN * ph, v = WIN * ph * (1 - ph);
    const D = wins.reduce((s, x) => s + (x - e) ** 2, 0) / v;
    return { ok: true, k, tot, ph, e, D, df: k - 1, p: chi2Upper(D, k - 1) };
  }

  // ---------- state ----------
  let state = JSON.parse(JSON.stringify(DEFAULTS));
  let step = 6, model = null, playing = null;

  function readInt(el){ const s = el.value.trim(); return /^\d+$/.test(s) ? parseInt(s, 10) : NaN; }

  function validate(){
    const pm = [], lm = []; let ok = true;
    const years = [['b1', 'd1', 'w1', 'Year 1'], ['b2', 'd2', 'w2', 'Year 2']];
    for (const [bk, dk, wk, name] of years){
      const bEl = $(bk), dEl = $(dk); bEl.classList.remove('bad', 'warn'); dEl.classList.remove('bad', 'warn');
      const b = readInt(bEl), d = readInt(dEl), w = state.mode === 'claim' ? 0 : state[wk];
      if (isNaN(b)){ bEl.classList.add('bad'); pm.push(['err', `${name}: breaches must be a whole number ≥ 0.`]); ok = false; }
      if (isNaN(d)){ dEl.classList.add('bad'); pm.push(['err', `${name}: trading days must be a whole number ≥ 0.`]); ok = false; }
      if (!isNaN(b) && !isNaN(d)){
        if (b > d){ bEl.classList.add('bad'); pm.push(['err', `${name}: breaches (${b}) can't exceed trading days (${d}).`]); ok = false; }
        else if (d < 1 && w > 0){ dEl.classList.add('bad'); pm.push(['err', `${name}: needs at least 1 trading day unless its weight is 0.`]); ok = false; }
        if (d > 260){ dEl.classList.add('warn'); pm.push(['warn', `${name}: ${d} trading days is more than a calendar year of trading.`]); }
        state[bk] = b; state[dk] = d;
      }
    }
    if (state.mode !== 'claim' && state.w1 === 0 && state.w2 === 0) pm.push(['warn', 'Both weights are 0. The prior is reduced to the minimum Beta(0.5, 0.5).']);
    document.querySelectorAll('#logBody tr').forEach((tr, i) => {
      const inp = tr.querySelector('.bx'); inp.classList.remove('bad', 'warn');
      const v = readInt(inp), lab = `Days ${i * WIN + 1}–${i * WIN + WIN}`;
      if (isNaN(v)){ inp.classList.add('bad'); lm.push(['err', `${lab}: breaches must be a whole number ≥ 0.`]); ok = false; return; }
      if (v > WIN){ inp.classList.add('bad'); lm.push(['err', `${lab}: breaches (${v}) can't exceed the ${WIN} trading days in the window.`]); ok = false; return; }
      if (v > 5){ inp.classList.add('warn'); lm.push(['warn', `${lab}: ${v} breaches in 25 days is extreme for a 99% VaR model.`]); }
      state.wins[i] = v;
    });
    renderMsgs($('priorMsgs'), pm); renderMsgs($('logMsgs'), lm);
    return ok;
  }
  function renderMsgs(el, list){ el.innerHTML = list.map(([t, m]) => `<div class="msg ${t}">${m}</div>`).join(''); }

  function compute(){
    const { b1, d1, b2, d2, w1, w2, mode, claimN } = state;
    let a0raw, b0raw;
    if (mode === 'claim'){ a0raw = claimN * CLAIM; b0raw = claimN * (1 - CLAIM); }
    else { a0raw = w1 * b1 + w2 * b2; b0raw = w1 * (d1 - b1) + w2 * (d2 - b2); }
    const a0 = Math.max(0.5, a0raw), b0 = Math.max(0.5, b0raw);
    const states = [];
    let B = 0, D = 0;
    for (let k = 0; k <= NW; k++){
      let wb = null;
      if (k > 0){ wb = state.wins[k - 1]; B += wb; D += WIN; }
      const a = a0 + B, b = b0 + (D - B), p = a / (a + b), R = DAYS - D;
      const varTheta = p * (1 - p) / (a + b + 1);
      const luck = R * p * (1 - p), ign = R * R * varTheta;
      const pmf = R > 0 ? betaBinomPmf(R, a, b) : [1];
      const pts = pmf.map((q, y) => [B + y, q]);
      let pg = 0, py = 0, pr = 0, mean = 0, emult = 0, c = 0, lo = null, hi = null;
      pts.forEach(([n, q]) => {
        const z = zoneOf(n); if (z === 'g') pg += q; else if (z === 'y') py += q; else pr += q;
        mean += n * q; emult += MULT(n) * q; c += q;
        if (lo === null && c >= 0.025) lo = n; if (hi === null && c >= 0.975) hi = n;
      });
      if (hi === null) hi = B + R;
      const ci = [betaQuantile(0.025, a, b), betaQuantile(0.975, a, b)];
      const pAbove = 1 - betaCdf(CLAIM, a, b);
      const disp = dispersion(state.wins.slice(0, k));
      states.push({ k, D, wb, B, a, b, p, R, pmf, pts, pg, py, pr, mean, emult, lo, hi, ci, pAbove, luck, ign, varTheta, disp,
        winRate: k ? wb / WIN : null, ytd: D ? B / D : null, zoneNow: zoneOf(B) });
    }
    return { a0, b0, a0raw, b0raw, priorMean: a0 / (a0 + b0), priorN: a0 + b0, states };
  }

  // ---------- rendering ----------
  function renderPrior(m){
    $('r1').textContent = state.d1 ? `${pc(state.b1 / state.d1)} breach rate` : '— breach rate';
    $('r2').textContent = state.d2 ? `${pc(state.b2 / state.d2)} breach rate` : '— breach rate';
    $('w1v').textContent = f2(state.w1); $('w2v').textContent = f2(state.w2); $('cnv').textContent = state.claimN;
    const claim = state.mode === 'claim';
    $('weightSliders').hidden = claim; $('claimSlider').hidden = !claim; $('claimNote').hidden = !claim;
    $('yearGrid').classList.toggle('dim', claim);
    $('weightSummary').textContent = claim
      ? `Model's claim 1% · ${state.claimN} pseudo-days`
      : `Y1 ${f2(state.w1)} · Y2 ${f2(state.w2)} → ${Math.round(m.priorN)} pseudo-days`;
    $('priorResult').textContent = `Prior θ ~ Beta(${f1(m.a0)}, ${f1(m.b0)}) · mean ${pc(m.priorMean)} · worth ${Math.round(m.priorN)} trading days`;
    document.querySelectorAll('.chip').forEach(c => {
      const p = c.dataset.preset;
      const on = claim ? p === 'claim' :
        (p === 'weighted' && state.w1 === 0.4 && state.w2 === 0.5) || (p === 'pooled' && state.w1 === 1 && state.w2 === 1) || (p === 'recent' && state.w1 === 0 && state.w2 === 1);
      c.classList.toggle('on', on);
    });
  }

  const zoneChip = z => `<span class="zchip ${z}">${ZNAME[z]}</span>`;

  function renderKpis(m){
    const s = m.states[step];
    $('dayLabel').textContent = step === 0 ? 'Preseason (prior only)' : `After day ${s.D}`;
    if (step === 0){
      $('k1').textContent = '—'; $('k1s').textContent = 'no days yet';
      $('k2').textContent = '—'; $('k2s').textContent = '0 / 0 days';
    } else {
      $('k1').textContent = pc(s.winRate, 1); $('k1s').textContent = `${s.wb} / ${WIN} days`;
      $('k2').textContent = pc(s.ytd); $('k2s').innerHTML = `${s.B} / ${s.D} days · now ${zoneChip(s.zoneNow)}`;
    }
    $('k3bar').innerHTML = `<div class="g" style="width:${s.pg * 100}%"></div><div class="y" style="width:${s.py * 100}%"></div><div class="r" style="width:${s.pr * 100}%"></div>`;
    $('k3').innerHTML = `<span class="g"><small>GREEN</small>${pc0(s.pg)}</span><span class="y"><small>YELLOW</small>${pc0(s.py)}</span><span class="r"><small>RED</small>${pc0(s.pr)}</span>`;
    $('k3s').textContent = s.R ? `Expected ${f1(s.mean)} breaches · 95%: ${s.lo}–${s.hi}` : `Year complete: ${s.B} breaches, ${ZNAME[s.zoneNow].toLowerCase()} zone`;
    $('k4').textContent = pc(s.p);
    $('k4s').textContent = `P(θ > 1%) = ${pc0(s.pAbove)}`;
    $('stepRange').value = step;
  }

  function renderCluster(m){
    const s = m.states[step], d = s.disp, fl = $('clusterFlag'), tx = $('clusterTxt');
    if (!d.ok){
      fl.className = 'flag mono na'; fl.textContent = 'Not enough data';
      tx.textContent = `The check needs at least 3 windows and 3 breaches. So far: ${d.k} window${d.k === 1 ? '' : 's'}, ${d.tot} breach${d.tot === 1 ? '' : 'es'}.`;
      return;
    }
    const cls = d.p < 0.05 ? 'hot' : d.p < 0.10 ? 'maybe' : 'ok';
    fl.className = `flag mono ${cls}`;
    fl.textContent = cls === 'hot' ? 'Clustering detected' : cls === 'maybe' ? 'Possible clustering' : 'No clear clustering';
    const msg = cls === 'ok'
      ? 'Breaches are spread about as evenly across windows as independent breaches would be.'
      : 'Breaches are bunched into a few windows more than independent breaches would be. The Beta-Binomial forecast assumes independence, so treat its range as too narrow.';
    tx.innerHTML = `Dispersion across ${d.k} windows: D = ${f1(d.D)} on ${d.df} degrees of freedom, p = ${d.p.toFixed(3)}. ${msg} See Guide, section 6.`;
  }

  function renderSteps(m){
    const s = m.states[step], ol = $('steps'), items = [];
    if (state.mode === 'claim'){
      items.push(`<li><span class="what">Build the prior from the model's 1% claim.</span>
<pre>α₀ = ${state.claimN} × 0.01 = ${f1(m.a0)}
β₀ = ${state.claimN} × 0.99 = ${f1(m.b0)}
prior mean = <span class="res">${pc(m.priorMean)}</span></pre></li>`);
    } else {
      items.push(`<li><span class="what">Build the prior from two backtest years.</span>
<pre>α₀ = ${f2(state.w1)}×${state.b1} + ${f2(state.w2)}×${state.b2} = ${f1(m.a0)}
β₀ = ${f2(state.w1)}×${state.d1 - state.b1} + ${f2(state.w2)}×${state.d2 - state.b2} = ${f1(m.b0)}
prior mean = ${f1(m.a0)} / ${f1(m.a0 + m.b0)} = <span class="res">${pc(m.priorMean)}</span></pre></li>`);
    }
    if (step > 0){
      items.push(`<li><span class="what">Update with every day so far: add breaches to α and clean days to β.</span>
<pre>α = ${f1(m.a0)} + ${s.B} breaches = ${f1(s.a)}
β = ${f1(m.b0)} + ${s.D - s.B} clean days = ${f1(s.b)}
θ̂ = ${f1(s.a)} / ${f1(s.a + s.b)} = <span class="res">${pc(s.p)}</span>   P(θ > 1%) = ${pc0(s.pAbove)}</pre></li>`);
    } else {
      items.push(`<li><span class="what">No days yet, so the posterior equals the prior.</span>
<pre>θ̂ = <span class="res">${pc(s.p)}</span>   P(θ > 1%) = ${pc0(s.pAbove)}</pre></li>`);
    }
    if (s.R > 0){
      items.push(`<li><span class="what">Item 3: forecast the year-end zone. Remaining breaches ~ Beta-Binomial(${s.R}, ${f1(s.a)}, ${f1(s.b)}).</span>
<pre>R = ${DAYS} − ${s.D} = ${s.R} days left
expected total = ${s.B} + ${s.R} × ${(s.p).toFixed(4)} = ${f1(s.mean)}
P(green ≤4) = ${pc0(s.pg)}  P(yellow 5–9) = ${pc0(s.py)}  P(red ≥10) = ${pc0(s.pr)}
expected multiplier = <span class="res">×${s.emult.toFixed(2)}</span></pre></li>`);
      const sd = Math.sqrt(s.luck + s.ign);
      items.push(`<li><span class="what">Item 4: the uncertainty curve and where its width comes from.</span>
<pre>luck      = R·p(1−p)  = ${s.luck.toFixed(2)}
ignorance = R²·Var(θ) = ${s.R}² × ${s.varTheta.toExponential(2)} = ${s.ign.toFixed(2)}
SD = √${(s.luck + s.ign).toFixed(2)} = ${sd.toFixed(2)} breaches
95% range (exact) = <span class="res">${s.lo} – ${s.hi} breaches</span></pre></li>`);
    } else {
      items.push(`<li><span class="what">Year complete. R = 0, so there is nothing left to forecast.</span>
<pre>final count = <span class="res">${s.B} breaches → ${ZNAME[s.zoneNow]} (×${MULT(s.B).toFixed(2)})</span>
θ̂ = ${pc(s.p)} still carries uncertainty: 95% credible range ${pc(s.ci[0])}–${pc(s.ci[1])}.</pre></li>`);
    }
    ol.innerHTML = items.join('');
  }

  const ns = 'http://www.w3.org/2000/svg';
  const el = (tag, attrs, txt) => { const e = document.createElementNS(ns, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); if (txt != null) e.textContent = txt; return e; };
  const ZC = { g: '#6fd3a0', y: '#f2c46d', r: '#f2818a' };

  function curveMax(m){
    let hi = 12; m.states.forEach(s => { hi = Math.max(hi, s.hi + 2, s.B + 2); });
    return Math.min(DAYS, hi);
  }

  function renderCurve(m){
    const svg = $('curveChart'); svg.innerHTML = '';
    const W = 640, Hh = 300, L = 40, Rr = 14, T = 30, Bm = 34;
    const xmax = curveMax(m), slot = (W - L - Rr) / (xmax + 1);
    const X = n => L + (n + 0.5) * slot;
    const cur = m.states[step];
    let ymax = 0; cur.pts.forEach(([n, q]) => { if (n <= xmax) ymax = Math.max(ymax, q); });
    for (let k = 0; k < step; k++) m.states[k].pts.forEach(([n, q]) => { if (n <= xmax) ymax = Math.max(ymax, q); });
    ymax = Math.min(1, Math.max(0.05, ymax * 1.15));
    const Y = v => Hh - Bm - v / ymax * (Hh - Bm - T);
    // zone backgrounds
    [['g', 0, 4], ['y', 5, 9], ['r', 10, xmax]].forEach(([z, a, b]) => {
      if (a > xmax) return;
      const x0 = L + a * slot, x1 = L + (Math.min(b, xmax) + 1) * slot;
      svg.appendChild(el('rect', { x: x0, y: T - 8, width: x1 - x0, height: Hh - Bm - T + 8, fill: ZC[z], opacity: 0.05 }));
      svg.appendChild(el('text', { x: (x0 + x1) / 2, y: T - 14, 'text-anchor': 'middle', style: `fill:${ZC[z]};font-weight:600;letter-spacing:.06em` }, ZNAME[z].toUpperCase()));
    });
    // y grid
    const yt = ymax > 0.5 ? 0.2 : ymax > 0.2 ? 0.1 : 0.05;
    for (let v = 0; v <= ymax + 1e-9; v += yt){
      svg.appendChild(el('line', { x1: L, x2: W - Rr, y1: Y(v), y2: Y(v), stroke: '#1c2638' }));
      svg.appendChild(el('text', { x: L - 6, y: Y(v) + 4, 'text-anchor': 'end' }, Math.round(v * 100) + '%'));
    }
    svg.appendChild(el('line', { x1: L, x2: W - Rr, y1: Hh - Bm, y2: Hh - Bm, stroke: '#2a3850' }));
    const every = xmax > 30 ? 5 : xmax > 16 ? 2 : 1;
    for (let n = 0; n <= xmax; n += every) svg.appendChild(el('text', { x: X(n), y: Hh - 16, 'text-anchor': 'middle' }, n));
    svg.appendChild(el('text', { x: (L + W - Rr) / 2, y: Hh - 2, 'text-anchor': 'middle' }, 'breaches at day 250'));
    // bars for current window
    const bw = Math.max(2, slot * 0.72);
    cur.pts.forEach(([n, q]) => {
      if (n > xmax || q < 1e-5) return;
      const inBand = n >= cur.lo && n <= cur.hi;
      svg.appendChild(el('rect', { x: X(n) - bw / 2, y: Y(q), width: bw, height: Hh - Bm - Y(q), rx: 2, fill: ZC[zoneOf(n)], opacity: inBand ? 0.9 : 0.4 }));
    });
    const line = pts => pts.filter(([n]) => n <= xmax).map(([n, q], i) => `${i ? 'L' : 'M'}${X(n).toFixed(1)},${Math.max(T - 4, Y(q)).toFixed(1)}`).join('');
    for (let k = 1; k < step; k++) svg.appendChild(el('path', { d: line(m.states[k].pts), fill: 'none', stroke: 'rgba(230,235,242,.22)', 'stroke-width': 1.2 }));
    if (step > 0) svg.appendChild(el('path', { d: line(m.states[0].pts), fill: 'none', stroke: '#8a96a8', 'stroke-width': 1.5, 'stroke-dasharray': '5 4' }));
    // breaches so far
    if (step > 0 && cur.R > 0){
      const bx = L + cur.B * slot;
      svg.appendChild(el('line', { x1: bx, x2: bx, y1: T + 4, y2: Hh - Bm, stroke: '#7fb2ff', 'stroke-width': 1.5, 'stroke-dasharray': '1 3' }));
      svg.appendChild(el('text', { x: bx - 4, y: T + 14, 'text-anchor': 'end', style: 'fill:#7fb2ff' }, `${cur.B} so far`));
    }
    const label = cur.R > 0 ? `E = ${f1(cur.mean)} · 95%: ${cur.lo}–${cur.hi}` : `Final: ${cur.B} breaches (${ZNAME[cur.zoneNow]})`;
    svg.appendChild(el('text', { x: W - Rr, y: T + 10, 'text-anchor': 'end', style: 'fill:#e6ebf2;font-weight:600' }, label));
    // variance split
    const vs = $('varSplit');
    if (cur.R > 0){
      const tot = cur.luck + cur.ign, lp = cur.luck / tot * 100;
      vs.innerHTML = `<div class="muted">Where the uncertainty comes from (variance in remaining breaches)</div>
        <div class="var-bar"><div class="luck" style="width:${lp}%"></div><div class="ign" style="width:${100 - lp}%"></div></div>
        <div class="var-labels"><span><span style="color:var(--blue)">■</span> Luck <b>${cur.luck.toFixed(2)}</b> (${lp.toFixed(0)}%)</span><span><span style="color:#c9a3ff">■</span> Rate uncertainty <b>${cur.ign.toFixed(2)}</b> (${(100 - lp).toFixed(0)}%)</span><span>SD <b>${Math.sqrt(tot).toFixed(2)}</b> breaches</span></div>`;
    } else vs.innerHTML = `<div class="muted">All 250 days are in. The year-end count is known, so luck and rate uncertainty about it are both zero.</div>`;
  }

  function renderZones(m){
    const svg = $('zoneChart'); svg.innerHTML = '';
    const W = 640, Hh = 200, L = 40, Rr = 14, T = 10, Bm = 28;
    const X = d => L + d / DAYS * (W - L - Rr), Y = v => T + (1 - v) * (Hh - T - Bm);
    const S = m.states;
    const band = (lo, hi, col) => {
      const up = S.map(s => `${X(s.D).toFixed(1)},${Y(hi(s)).toFixed(1)}`), dn = S.slice().reverse().map(s => `${X(s.D).toFixed(1)},${Y(lo(s)).toFixed(1)}`);
      svg.appendChild(el('path', { d: `M${up.join('L')}L${dn.join('L')}Z`, fill: col, opacity: 0.85 }));
    };
    band(s => 0, s => s.pg, ZC.g);
    band(s => s.pg, s => s.pg + s.py, ZC.y);
    band(s => s.pg + s.py, s => 1, ZC.r);
    [0, 0.5, 1].forEach(v => svg.appendChild(el('text', { x: L - 6, y: Y(v) + 4, 'text-anchor': 'end' }, Math.round(v * 100) + '%')));
    svg.appendChild(el('line', { x1: L, x2: W - Rr, y1: Y(0.5), y2: Y(0.5), stroke: '#0a0e16', 'stroke-dasharray': '3 3', opacity: .6 }));
    [0, 50, 100, 150, 200, 250].forEach(d => svg.appendChild(el('text', { x: X(d), y: Hh - 8, 'text-anchor': 'middle' }, d === 0 ? 'Pre' : `D${d}`)));
    if (step < NW) svg.appendChild(el('rect', { x: X(S[step].D), y: T, width: X(DAYS) - X(S[step].D), height: Hh - T - Bm, fill: 'rgba(10,14,22,.62)' }));
    const cx = X(S[step].D);
    svg.appendChild(el('line', { x1: cx, x2: cx, y1: T, y2: Hh - Bm, stroke: '#e6ebf2', 'stroke-width': 1.5 }));
    addClick(svg, W, L, Rr, T, Hh - T - Bm);
  }

  function renderTraj(m){
    const svg = $('trajChart'); svg.innerHTML = '';
    const W = 640, Hh = 260, L = 46, Rr = 14, T = 14, Bm = 30;
    const S = m.states;
    let ymax = 0.03; S.forEach(s => { ymax = Math.max(ymax, s.ci[1], s.ytd || 0); });
    ymax = Math.ceil(ymax * 1.1 * 200) / 200;
    const X = d => L + d / DAYS * (W - L - Rr), Y = v => T + (1 - v / ymax) * (Hh - T - Bm);
    const st = ymax > 0.06 ? 0.02 : 0.01;
    for (let v = 0; v <= ymax + 1e-9; v += st){
      svg.appendChild(el('line', { x1: L, x2: W - Rr, y1: Y(v), y2: Y(v), stroke: '#1c2638' }));
      svg.appendChild(el('text', { x: L - 6, y: Y(v) + 4, 'text-anchor': 'end' }, Math.round(v * 100) + '%'));
    }
    [0, 50, 100, 150, 200, 250].forEach(d => svg.appendChild(el('text', { x: X(d), y: Hh - 10, 'text-anchor': 'middle' }, d === 0 ? 'Pre' : `D${d}`)));
    const up = S.map(s => `${X(s.D).toFixed(1)},${Y(s.ci[1]).toFixed(1)}`), dn = S.slice().reverse().map(s => `${X(s.D).toFixed(1)},${Y(s.ci[0]).toFixed(1)}`);
    svg.appendChild(el('path', { d: `M${up.join('L')}L${dn.join('L')}Z`, fill: 'rgba(111,211,160,.13)' }));
    svg.appendChild(el('line', { x1: L, x2: W - Rr, y1: Y(CLAIM), y2: Y(CLAIM), stroke: '#e6ebf2', 'stroke-dasharray': '6 4', opacity: .55 }));
    svg.appendChild(el('text', { x: L + 6, y: Y(CLAIM) - 5, 'text-anchor': 'start' }, '1% claim'));
    const ys = S.filter(s => s.ytd != null);
    svg.appendChild(el('path', { d: 'M' + ys.map(s => `${X(s.D).toFixed(1)},${Y(s.ytd).toFixed(1)}`).join('L'), fill: 'none', stroke: '#f2c46d', 'stroke-width': 1.8 }));
    svg.appendChild(el('path', { d: 'M' + S.map(s => `${X(s.D).toFixed(1)},${Y(s.p).toFixed(1)}`).join('L'), fill: 'none', stroke: '#6fd3a0', 'stroke-width': 2.6 }));
    S.forEach(s => {
      if (s.winRate == null) return;
      const cx = X(s.D - WIN / 2), op = s.k <= step ? 0.95 : 0.3;
      if (s.winRate > ymax){
        svg.appendChild(el('path', { d: `M${cx - 5},${T + 8}L${cx + 5},${T + 8}L${cx},${T}Z`, fill: '#7fb2ff', opacity: op }));
        svg.appendChild(el('text', { x: cx, y: T + 20, 'text-anchor': 'middle', style: 'fill:#7fb2ff', opacity: op }, pc(s.winRate, 0)));
      } else svg.appendChild(el('circle', { cx, cy: Y(s.winRate), r: 3.4, fill: '#7fb2ff', opacity: op }));
    });
    if (step < NW) svg.appendChild(el('rect', { x: X(S[step].D), y: T, width: X(DAYS) - X(S[step].D), height: Hh - T - Bm, fill: 'rgba(10,14,22,.55)' }));
    const cx = X(S[step].D);
    svg.appendChild(el('line', { x1: cx, x2: cx, y1: T, y2: Hh - Bm, stroke: '#e6ebf2', 'stroke-width': 1, opacity: .5 }));
    svg.appendChild(el('circle', { cx, cy: Y(S[step].p), r: 5, fill: '#6fd3a0', stroke: '#0a0e16', 'stroke-width': 2 }));
    addClick(svg, W, L, Rr, T, Hh - T - Bm);
  }

  function addClick(svg, W, L, Rr, T, h){
    const hit = el('rect', { x: L, y: T, width: W - L - Rr, height: h, fill: 'transparent', style: 'cursor:pointer' });
    hit.addEventListener('click', ev => {
      const r = svg.getBoundingClientRect(), px = (ev.clientX - r.left) / r.width * W;
      setStep(Math.round(Math.min(DAYS, Math.max(0, (px - L) / (W - L - Rr) * DAYS)) / WIN));
    });
    svg.appendChild(hit);
  }

  function buildTable(){
    const tb = $('logBody'); tb.innerHTML = '';
    state.wins.forEach((v, i) => {
      const tr = document.createElement('tr'); tr.dataset.k = i + 1;
      tr.innerHTML = `<td>${i * WIN + 1}–${i * WIN + WIN}</td>
        <td><input class="num bx" inputmode="numeric" value="${v}" aria-label="Breaches days ${i * WIN + 1} to ${i * WIN + WIN}"></td>
        <td class="c-win"></td><td class="c-ytd"></td><td class="c-th bayes"></td><td class="c-z"></td><td class="c-e"></td>`;
      tr.addEventListener('click', e => { if (e.target.tagName !== 'INPUT') setStep(i + 1); });
      tb.appendChild(tr);
    });
    tb.querySelectorAll('input').forEach(inp => inp.addEventListener('input', update));
  }
  function renderTable(m){
    document.querySelectorAll('#logBody tr').forEach((tr, i) => {
      const s = m.states[i + 1];
      tr.querySelector('.c-win').textContent = pc(s.winRate, 0);
      tr.querySelector('.c-ytd').textContent = pc(s.ytd);
      tr.querySelector('.c-th').textContent = pc(s.p);
      tr.querySelector('.c-z').innerHTML = s.R
        ? `<span class="mini-zone" title="G ${pc0(s.pg)} · Y ${pc0(s.py)} · R ${pc0(s.pr)}"><i class="g" style="width:${s.pg * 100}%"></i><i class="y" style="width:${s.py * 100}%"></i><i class="r" style="width:${s.pr * 100}%"></i></span>`
        : zoneChip(s.zoneNow);
      tr.querySelector('.c-e').textContent = s.R ? `${f1(s.mean)} (${s.lo}–${s.hi})` : `${s.B} final`;
      tr.classList.toggle('sel', i + 1 === step);
      tr.classList.toggle('future', i + 1 > step);
    });
  }

  function renderGuide(m){
    const ex = step >= 1 && step <= NW - 1 ? step : (step === 0 ? 1 : NW - 1);
    const s = m.states[ex], s0 = m.states[0], last = m.states[NW];
    const { b1, d1, b2, d2, w1, w2 } = state;
    const claim = state.mode === 'claim';
    const note = ex === step
      ? `Worked numbers below use your current inputs and the window selected in the Simulator (day ${s.D}).`
      : step === 0
        ? `Preseason is selected in the Simulator. Worked numbers from section 2 onward use day ${s.D}, the first window with data. Select a window from day 25 to 225 to change it.`
        : `Worked numbers below use your current inputs at day ${s.D}. Select a window from day 25 to 225 in the Simulator to change it.`;
    const priorDesc = claim
      ? `a prior centered on the model's own 1% claim, worth ${state.claimN} trading days`
      : `his two previous backtest years, weighted ${f2(w1)} and ${f2(w2)}, which put the breach rate at <b>${pc(m.priorMean)}</b>`;
    const q1 = step === 0
      ? `It's the start of the year, so the desk has no current-year days yet and there's nothing to update on. The Bayes forecast rests entirely on the prior: ${priorDesc}. That alone gives a <b>${pc0(s0.pg)}</b> chance of finishing the year green, ${pc0(s0.py)} yellow, and ${pc0(s0.pr)} red, with no posterior yet.`
      : `After ${s.D} trading days, the 99% VaR model has been breached <b>${s.B}</b> time${s.B === 1 ? '' : 's'} (${pc(s.ytd)} of days), so the traffic light currently reads ${ZNAME[s.zoneNow].toLowerCase()}. The prior comes from ${priorDesc}.`;
    const wShare = s.D / (s.D + m.priorN);
    const zoneNowTxt = s.zoneNow === 'g' && (s.py + s.pr) > 0.5
      ? ` The traffic light still reads green, yet Bayes puts a ${pc0(s.py + s.pr)} chance on finishing yellow or red. Counting only looks back, while the forecast looks ahead.`
      : '';
    const disp = last.disp;
    const dispTxt = disp.ok
      ? `For the full year in your inputs: ${disp.k} windows, ${disp.tot} breaches, expected ${disp.e.toFixed(2)} per window, D = ${f1(disp.D)} on ${disp.df} degrees of freedom, p = ${disp.p.toFixed(3)}. ${disp.p < 0.05 ? 'That is strong evidence of clustering.' : disp.p < 0.10 ? 'That is weak evidence of clustering.' : 'That shows no clear clustering.'}`
      : 'Your full-year inputs have fewer than 3 breaches, so the check has nothing to work with.';
    const tot = s.luck + s.ign, sd = Math.sqrt(tot);

    $('guideDyn').innerHTML = `
      <section class="card">
        <h2 class="g">1. The Question</h2>
        <p class="guide-note">${note}</p>
        <p>${q1} Where will the breach count end up at day 250, which Basel zone does that put the desk in, and how sure can we be?</p>
        <p>Bayes' theorem answers with one rule: <b>posterior ∝ likelihood × prior</b>. We start with what we believed about the model (the prior), weigh it against this year's breaches (the likelihood), and end up with an updated belief (the posterior). That posterior becomes the prior for the next 25 days.</p>
        <div class="callout"><b>The traffic light counts; Bayes forecasts.</b> The Basel test looks at the 250-day count only once the year is in. A Bayesian forecast gives the probability of each zone at any point in the year, so validators see trouble coming before the count crosses a line.</div>
      </section>

      <section class="card">
        <h2 class="g">2. Building the Prior</h2>
        <p>The model's true chance of a breach on any day, θ, is unknown. A correct 99% VaR model has θ = 1%. We describe our belief about θ with a <b>Beta distribution</b>, whose two parameters read as pseudo-breaches and pseudo-clean days:</p>
        <pre class="eq">α₀ = w₁·B₁ + w₂·B₂        β₀ = w₁·(D₁−B₁) + w₂·(D₂−B₂)</pre>
        <p>With your inputs (Year 1: ${b1} breaches in ${d1} days, Year 2: ${b2} in ${d2}, weights ${f2(w1)} and ${f2(w2)}), the history-based prior is Beta(${f1(Math.max(0.5, w1 * b1 + w2 * b2))}, ${f1(Math.max(0.5, w1 * (d1 - b1) + w2 * (d2 - b2)))}).${claim ? ` You've selected the model's-claim prior instead: α₀ = ${state.claimN} × 0.01 and β₀ = ${state.claimN} × 0.99, giving Beta(${f1(m.a0)}, ${f1(m.b0)}).` : ''}</p>
        <h3 class="g3">Two ways to set the prior</h3>
        <ul>
          <li><b>Track record (default).</b> Use the model's own last two years. It's honest about how the model has actually behaved, and recent years count more because markets and portfolios change.</li>
          <li><b>Model's claim.</b> Center the prior on exactly 1% and choose how many days of trust it deserves. This is the validator's stance: assume the model works until the data says otherwise.</li>
        </ul>
        <div class="callout"><b>Why the prior is so light.</b> A correct model produces only about 2.5 breaches a year. Two years of history hold just ${b1 + b2} breaches, so even a well-documented prior is weak, and each new breach moves the estimate noticeably. Kupiec (1995) made the same point: tail probabilities are inherently hard to verify.</div>
      </section>

      <section class="card">
        <h2 class="g">3. Updating Every 25 Days</h2>
        <p>Each window adds its breaches to α and its clean days to β. Because the Beta prior is conjugate to the Binomial, the posterior is again a Beta:</p>
        <pre class="eq">θ | data ~ Beta(α₀ + B, β₀ + D − B)
θ̂ = (α₀ + B) / (α₀ + β₀ + D)</pre>
        <p>At day ${s.D}: α = ${f1(m.a0)} + ${s.B} = ${f1(s.a)}, β = ${f1(m.b0)} + ${s.D - s.B} = ${f1(s.b)}, so θ̂ = <b>${pc(s.p)}</b>, with a 95% credible range of ${pc(s.ci[0])}–${pc(s.ci[1])}. The probability that the true breach rate is above the 1% claim is <b>${pc0(s.pAbove)}</b>.</p>
        <p>The data now carries ${pc0(wShare)} of the weight (${s.D} real days versus ${Math.round(m.priorN)} pseudo-days of prior).</p>
        <div class="callout"><b>Yesterday's posterior is today's prior.</b> Updating once with all 250 days gives exactly the same answer as updating 10 times, 25 days at a time. The order of the breaches doesn't change the posterior, which is precisely why clustering (section 6) needs a separate check.</div>
      </section>

      <section class="card">
        <h2 class="g">4. Item 3: Forecasting the Year-End Zone</h2>
        <p>The breaches still to come over the remaining R = ${DAYS} − ${s.D} = ${s.R} days follow a Beta-Binomial distribution. Add them to the ${s.B} already counted to get the year-end total, then sort each possible total into a zone:</p>
        <table class="zt">
          <tr><th>Zone</th><th>Breaches in 250 days</th><th>Multiplier (1996)</th></tr>
          <tr><td>${zoneChip('g')}</td><td>0–4</td><td>3.00</td></tr>
          <tr><td>${zoneChip('y')}</td><td>5, 6, 7, 8, 9</td><td>3.40, 3.50, 3.65, 3.75, 3.85</td></tr>
          <tr><td>${zoneChip('r')}</td><td>10 or more</td><td>4.00</td></tr>
        </table>
        <p>At day ${s.D}, the expected total is ${s.B} + ${s.R} × ${(s.p).toFixed(4)} = <b>${f1(s.mean)}</b>, with a 95% range of ${s.lo}–${s.hi}. Summing the probabilities by zone gives <b>${pc0(s.pg)} green, ${pc0(s.py)} yellow, and ${pc0(s.pr)} red</b>, and an expected multiplier of ×${s.emult.toFixed(2)}.${zoneNowTxt}</p>
        <div class="callout"><b>Yellow isn't proof of failure.</b> The zone boundaries are set by a Binomial with θ = 1%: yellow starts where a correct model would show that many breaches or fewer at least 95% of the time, and red where it would do so 99.99% of the time. So even a perfect model lands in yellow roughly one year in ten. The Bayesian question is different: given this year's data, how likely is it that θ itself is above 1%?</div>
      </section>

      <section class="card">
        <h2 class="g">5. Item 4: The Uncertainty Curve</h2>
        <p>Remaining breaches Y follow a <b>Beta-Binomial</b> distribution, which is a Binomial whose breach probability is itself uncertain:</p>
        <pre class="eq">P(Y = y) = C(R, y) · B(y + α, R − y + β) / B(α, β),    y = 0 … R
Year-end total = B + Y</pre>
        <p>The app computes this probability for every possible y, using log-gamma functions to avoid overflow, and draws it as bars colored by zone. The 95% range comes from reading the 2.5% and 97.5% points off the cumulative distribution.</p>
        <pre class="eq">Var(Y) = R·p(1−p)           ← luck: days still to be traded
       + R²·Var(θ)          ← rate uncertainty: we don't know θ
where p = α/(α+β),  Var(θ) = p(1−p)/(α+β+1)</pre>
        <p>At day ${s.D}: luck = ${s.R} × ${(s.p).toFixed(4)} × ${(1 - s.p).toFixed(4)} = ${s.luck.toFixed(2)}; rate uncertainty = ${s.R}² × ${s.varTheta.toExponential(2)} = ${s.ign.toFixed(2)}. The total is ${tot.toFixed(2)}, so the SD is ${sd.toFixed(2)} breaches, and the exact 95% range is <b>${s.lo}–${s.hi}</b>.</p>
        <div class="callout"><b>Why rate uncertainty matters more here than in baseball.</b> With so few breaches, the posterior on θ stays wide all year. Early on, rate uncertainty can rival luck as a source of spread, which is the Bayesian version of Kupiec's warning that one year of data can barely tell a 1% model from a 2% model.</div>
      </section>

      <section class="card">
        <h2 class="g">6. The Clustering Caveat</h2>
        <p>Everything above assumes breaches are <b>independent</b>: a breach today says nothing about tomorrow. Real breaches don't behave that way. When volatility jumps, a VaR model calibrated to calmer data is breached several days in a row, then catches up. Breaches arrive in clusters.</p>
        <h3 class="g3">Why it matters</h3>
        <ul>
          <li><b>Counts hide timing.</b> Seven breaches spread across the year and seven breaches in three weeks give the same count, the same posterior, and the same Basel zone. Only the second points to a model that reacts too slowly to changing markets.</li>
          <li><b>The forecast range is too narrow.</b> Under clustering, breaches come in bunches, so the year-end count varies more than the Beta-Binomial assumes. Read the 95% range as a lower bound on the real uncertainty.</li>
        </ul>
        <h3 class="g3">The standard test</h3>
        <p>Christoffersen (1998) tests independence from the day-by-day sequence of breaches. It compares the chance of a breach on the day after a breach with the chance of a breach on the day after a clean day. If those differ, breaches are clustered. Combined with Kupiec's count test, it gives the conditional-coverage test that validators use. It needs daily data, which this app doesn't collect.</p>
        <h3 class="g3">What the app checks instead</h3>
        <p>The Clustering Check card runs a dispersion test on the 25-day window counts. If breaches were independent, each window's count would vary only as much as a Binomial(25, p̂) allows. Bunching shows up as extra spread:</p>
        <pre class="eq">D = Σ (xᵢ − 25·p̂)² / (25·p̂·(1−p̂))   ~ χ²(k − 1) under independence</pre>
        <p>${dispTxt}</p>
        <div class="callout"><b>A blunt instrument.</b> With only a handful of breaches, this test has little power, and the χ² approximation is rough. A quiet flag doesn't prove independence. A loud one, though, is a strong hint to pull the daily data and run the Christoffersen test.</div>
      </section>

      <section class="card">
        <h2 class="g">7. What the Bayesian View Adds</h2>
        <ul>
          <li><b>A probability, not a verdict.</b> "There's a ${pc0(last.pAbove)} chance the true breach rate is above 1%" is easier to act on than "yellow zone."</li>
          <li><b>An early warning.</b> Zone probabilities update every 25 days, so a risk committee sees a likely yellow before day 250.</li>
          <li><b>Explicit trust.</b> The prior makes clear how much benefit of the doubt the model is getting, and the presets show how much that choice matters.</li>
        </ul>
      </section>`;
  }

  function renderAll(){
    if (!model) return;
    renderPrior(model); renderKpis(model); renderCluster(model); renderSteps(model);
    renderCurve(model); renderZones(model); renderTraj(model); renderTable(model); renderGuide(model);
  }
  function update(){ if (validate()) model = compute(); renderAll(); }
  function setStep(k){ step = Math.max(0, Math.min(NW, k)); renderAll(); }

  // ---------- wiring ----------
  function syncInputs(){
    ['b1', 'd1', 'b2', 'd2'].forEach(k => $(k).value = state[k]);
    $('w1').value = state.w1; $('w2').value = state.w2; $('cn').value = state.claimN;
  }
  ['b1', 'd1', 'b2', 'd2'].forEach(k => $(k).addEventListener('input', update));
  $('w1').addEventListener('input', e => { state.w1 = +e.target.value; state.mode = 'history'; update(); });
  $('w2').addEventListener('input', e => { state.w2 = +e.target.value; state.mode = 'history'; update(); });
  $('cn').addEventListener('input', e => { state.claimN = +e.target.value; update(); });
  document.querySelectorAll('.chip').forEach(c => c.addEventListener('click', () => {
    const p = c.dataset.preset;
    if (p === 'claim') state.mode = 'claim';
    else {
      state.mode = 'history';
      if (p === 'weighted'){ state.w1 = 0.4; state.w2 = 0.5; }
      if (p === 'pooled'){ state.w1 = 1; state.w2 = 1; }
      if (p === 'recent'){ state.w1 = 0; state.w2 = 1; }
    }
    syncInputs(); update();
  }));
  $('resetBtn').addEventListener('click', () => {
    stopPlay(); state = JSON.parse(JSON.stringify(DEFAULTS)); syncInputs(); buildTable(); step = 6; update();
  });
  $('stepRange').addEventListener('input', e => { stopPlay(); setStep(+e.target.value); });
  $('prevBtn').addEventListener('click', () => { stopPlay(); setStep(step - 1); });
  $('nextBtn').addEventListener('click', () => { stopPlay(); setStep(step + 1); });
  function stopPlay(){ if (playing){ clearInterval(playing); playing = null; $('playBtn').textContent = 'Play'; } }
  $('playBtn').addEventListener('click', () => {
    if (playing) return stopPlay();
    if (step >= NW) setStep(0);
    $('playBtn').textContent = 'Pause';
    playing = setInterval(() => { if (step >= NW) return stopPlay(); setStep(step + 1); }, 900);
  });
  document.querySelectorAll('.tab').forEach(t => t.addEventListener('click', () => {
    document.querySelectorAll('.tab').forEach(x => { x.classList.toggle('active', x === t); x.setAttribute('aria-selected', x === t); });
    $('tab-sim').classList.toggle('hidden', t.dataset.tab !== 'sim');
    $('tab-guide').classList.toggle('hidden', t.dataset.tab !== 'guide');
  }));

  buildTable(); update();
})();
