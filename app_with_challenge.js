/* Lead–Lag Design Lab
   Mirrors MATLAB:
   G = 1/(s*(s+1)*(s+2))
   C = K*(s+a)/s * N*(s+b)/(s+b*N)
   y = lsim(feedback(G*C,1), r-n) + lsim(feedback(G,C), l)
   u = lsim(feedback(C,G), r-n) + lsim(feedback(G*C,1), l)
*/


function polyEvalComplexDesc(coeffs, z) {
  // coeffs descending powers, complex z = {re,im}
  let y = c(0, 0);
  for (let i = 0; i < coeffs.length; i++) {
    y = cAdd(cMul(y, z), c(coeffs[i], 0));
  }
  return y;
}

function durandKernerRootsRealPoly(den, maxIter = 200, tol = 1e-10) {
  // den: real coeffs descending powers [a0..an], a0 != 0
  if (!Array.isArray(den) || den.length < 2) return [];
  const n = den.length - 1;
  let a0 = den[0];
  if (!Number.isFinite(a0) || Math.abs(a0) < 1e-14) return [];

  // Normalize to monic for better conditioning
  const p = den.map(v => v / a0);

  // Cauchy root bound: 1 + max |ai|
  let R = 1;
  for (let i = 1; i < p.length; i++) R = Math.max(R, 1 + Math.abs(p[i]));

  // Initial guesses on a circle
  let roots = [];
  for (let k = 0; k < n; k++) {
    const ang = 2 * Math.PI * k / n;
    roots.push(c(R * Math.cos(ang), R * Math.sin(ang)));
  }

  // Iterations
  for (let it = 0; it < maxIter; it++) {
    let maxDelta = 0;

    for (let i = 0; i < n; i++) {
      let zi = roots[i];
      let fzi = polyEvalComplexDesc(p, zi);

      // denom = Π_{j≠i} (zi - zj)
      let denom = c(1, 0);
      for (let j = 0; j < n; j++) {
        if (j === i) continue;
        denom = cMul(denom, cSub(zi, roots[j]));
      }

      // If denom is tiny, perturb slightly
      const denomAbs = cAbs(denom);
      if (denomAbs < 1e-14) {
        zi = c(zi.re + 1e-6, zi.im + 1e-6);
        roots[i] = zi;
        continue;
      }

      const delta = cDiv(fzi, denom);
      const ziNew = cSub(zi, delta);
      roots[i] = ziNew;
      maxDelta = Math.max(maxDelta, cAbs(delta));
    }

    if (maxDelta < tol) break;
  }

  // Clean tiny imaginary parts
  return roots.map(z => ({
    re: z.re,
    im: Math.abs(z.im) < 1e-8 ? 0 : z.im
  }));
}

function parseCoeffs(text) {
  const parts = text.split(/[\s,;]+/).map(s => s.trim()).filter(Boolean);
  const coeffs = parts.map(Number);
  if (coeffs.length === 0 || coeffs.some(x => !Number.isFinite(x))) {
    throw new Error("Invalid coefficient list.");
  }
  // remove leading zeros
  while (coeffs.length > 1 && Math.abs(coeffs[0]) < 1e-14) coeffs.shift();
  return coeffs;
}

function polyAdd(a, b) {
  const n = Math.max(a.length, b.length);
  const A = Array(n).fill(0);
  const B = Array(n).fill(0);

  for (let i = 0; i < a.length; i++) {
    A[n - a.length + i] = a[i];
  }
  for (let i = 0; i < b.length; i++) {
    B[n - b.length + i] = b[i];
  }

  const out = A.map((v, i) => v + B[i]);

  // --- REMOVE LEADING ZEROS ---
  while (out.length > 1 && Math.abs(out[0]) < 1e-14) {
    out.shift();
  }

  return out;
}
function polyMul(a, b) {
  const out = Array(a.length + b.length - 1).fill(0);
  for (let i = 0; i < a.length; i++) {
    for (let j = 0; j < b.length; j++) out[i + j] += a[i] * b[j];
  }
  // trim leading zeros
  while (out.length > 1 && Math.abs(out[0]) < 1e-14) out.shift();
  return out;
}

function polyScale(a, k) {
  return a.map(v => v * k);
}

// complex helpers as {re, im}
function c(re, im) { return { re, im }; }
function cAdd(z1, z2) { return c(z1.re + z2.re, z1.im + z2.im); }
function cSub(z1, z2) { return c(z1.re - z2.re, z1.im - z2.im); }
function cMul(z1, z2) { return c(z1.re * z2.re - z1.im * z2.im, z1.re * z2.im + z1.im * z2.re); }
function cDiv(z1, z2) {
  const d = z2.re * z2.re + z2.im * z2.im;
  return c((z1.re * z2.re + z1.im * z2.im) / d, (z1.im * z2.re - z1.re * z2.im) / d);
}
function cAbs(z) { return Math.hypot(z.re, z.im); }
function cArgDeg(z) { return Math.atan2(z.im, z.re) * 180 / Math.PI; }

function polyEvalComplex(coeffs, s) {
  // coeffs descending powers
  let z = c(0, 0);
  for (let i = 0; i < coeffs.length; i++) {
    z = cAdd(cMul(z, s), c(coeffs[i], 0));
  }
  return z;
}

function tfFreqResp(num, den, w) {
  const s = c(0, w); // j*w
  const n = polyEvalComplex(num, s);
  const d = polyEvalComplex(den, s);
  return cDiv(n, d);
}

function logspace(a, b, n) {
  const out = [];
  for (let i = 0; i < n; i++) {
    const t = i / (n - 1);
    out.push(Math.pow(10, a + (b - a) * t));
  }
  return out;
}

function unwrapPhaseDeg(ph) {
  // basic unwrap in degrees
  const out = ph.slice();
  for (let i = 1; i < out.length; i++) {
    let dp = out[i] - out[i - 1];
    while (dp > 180) { out[i] -= 360; dp -= 360; }
    while (dp < -180) { out[i] += 360; dp += 360; }
  }
  return out;
}

function interpLogW(x1, y1, x2, y2, yTarget) {
  // linear interpolation in log10(w) between (x1,y1) and (x2,y2)
  // returns w where y crosses yTarget (assuming y1 and y2 straddle)
  const lx1 = Math.log10(x1), lx2 = Math.log10(x2);
  const t = (yTarget - y1) / (y2 - y1);
  const lx = lx1 + t * (lx2 - lx1);
  return Math.pow(10, lx);
}

function computeMargins(w, mag, phaseDegUnwrapped) {
  // Gain crossover: |L| = 1  => find straddle in log space
  let Wcp = NaN, Pm = NaN;
  for (let i = 1; i < w.length; i++) {
    if ((mag[i - 1] - 1) * (mag[i] - 1) <= 0) {
      if (mag[i - 1] === mag[i]) continue;
      Wcp = interpLogW(w[i - 1], mag[i - 1], w[i], mag[i], 1);
      // interpolate phase at Wcp in log-w
      const t = (Math.log10(Wcp) - Math.log10(w[i - 1])) / (Math.log10(w[i]) - Math.log10(w[i - 1]));
      const ph = phaseDegUnwrapped[i - 1] + t * (phaseDegUnwrapped[i] - phaseDegUnwrapped[i - 1]);
      Pm = 180 + ph;
      break;
    }
  }

  // Phase crossover: phase = -180
  let Wcg = NaN, Gm = NaN;
  for (let i = 1; i < w.length; i++) {
    const p1 = phaseDegUnwrapped[i - 1], p2 = phaseDegUnwrapped[i];
    if ((p1 + 180) * (p2 + 180) <= 0) {
      if (p1 === p2) continue;
      Wcg = interpLogW(w[i - 1], p1, w[i], p2, -180);
      // interpolate magnitude at Wcg
      const t = (Math.log10(Wcg) - Math.log10(w[i - 1])) / (Math.log10(w[i]) - Math.log10(w[i - 1]));
      const m = mag[i - 1] + t * (mag[i] - mag[i - 1]);
      Gm = 1 / m;
      break;
    }
  }

  return { Gm, Pm, Wcg, Wcp };
}

// State-space from TF in controllable canonical form, simulate via RK4
function tfToSS(num, den) {
  // Normalize so den[0] = 1
  const a0 = den[0];
  const denN = den.map(v => v / a0);
  const numN = num.map(v => v / a0);

  const n = denN.length - 1; // order
  if (n <= 0) {
    // static gain
    const gain = (numN.length ? numN[numN.length - 1] : 0) / denN[denN.length - 1];
    return { A: [], B: [], C: [], D: gain };
  }

  // Pad numerator to length n+1 (descending powers)
  const numPad = Array(n + 1).fill(0);
  for (let i = 0; i < numN.length; i++) {
    numPad[n + 1 - numN.length + i] = numN[i];
  }

  // Companion A with correct coefficient ordering
  const A = Array.from({ length: n }, () => Array(n).fill(0));
  for (let i = 0; i < n - 1; i++) A[i][i + 1] = 1;

  // last row must be [-an, -a(n-1), ..., -a1]
  for (let j = 0; j < n; j++) A[n - 1][j] = -denN[n - j];

  const B = Array(n).fill(0);
  B[n - 1] = 1;

  // TF numerator: b0 s^n + b1 s^(n-1) + ... + bn
  const D = numPad[0];

  // C must match the same state ordering: [bn - an*D, ..., b1 - a1*D]
  const C = Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    const bi = numPad[n - i];     // bn, b(n-1), ..., b1
    const ai = denN[n - i];       // an, a(n-1), ..., a1
    C[i] = bi - ai * D;
  }

  return { A, B, C, D };
}


function rk4SimSS(ss, u, t) {
  const { A, B, C, D } = ss;
  const n = A.length;
  let x = Array(n).fill(0);
  const y = Array(t.length).fill(0);

  function AxPlusBu(xv, uv) {
    const out = Array(n).fill(0);
    for (let i = 0; i < n; i++) {
      let s = 0;
      for (let j = 0; j < n; j++) s += A[i][j] * xv[j];
      s += B[i] * uv;
      out[i] = s;
    }
    return out;
  }

  function addScaled(xv, kv, h) {
    const out = Array(n).fill(0);
    for (let i = 0; i < n; i++) out[i] = xv[i] + h * kv[i];
    return out;
  }

  for (let k = 0; k < t.length; k++) {
    const uk = u[k];
    // output
    let yk = D * uk;
    for (let i = 0; i < n; i++) yk += C[i] * x[i];
    y[k] = yk;

    if (k === t.length - 1) break;
    const h = t[k + 1] - t[k];

    const k1 = AxPlusBu(x, uk);
    const k2 = AxPlusBu(addScaled(x, k1, h / 2), u[k + 1] /* close enough */);
    const k3 = AxPlusBu(addScaled(x, k2, h / 2), u[k + 1]);
    const k4 = AxPlusBu(addScaled(x, k3, h), u[k + 1]);

    for (let i = 0; i < n; i++) {
      x[i] = x[i] + (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
  }
  return y;
}

function simulateTF(num, den, u, t) {
  const ss = tfToSS(num, den);
  return rk4SimSS(ss, u, t);
}

function closedLoopTF(numP, denP, numC, denC) {
  // P = numP/denP, C = numC/denC
  // L = P*C = (numP*numC)/(denP*denC)
  const numL = polyMul(numP, numC);
  const denL = polyMul(denP, denC);

  // T = L/(1+L) = numL / (denL + numL)
  const denCL = polyAdd(denL, numL);
  return { numL, denL, denCL };
}

function feedbackTF(numFwd, denFwd, numFb, denFb) {
  // feedback(Fwd, Fb) = Fwd / (1 + Fwd*Fb)
  const numLoop = polyMul(numFwd, numFb);
  const denLoop = polyMul(denFwd, denFb);
  const denOut = polyAdd(denLoop, numLoop);
  const numOut = polyMul(numFwd, denFb);
  const denOut2 = polyMul(denFwd, denOut); // denFwd * (denLoop+numLoop) / denFb? careful:
  // Derivation with polynomials:
  // Fwd = Nf/Df, Fb = Nb/Db
  // Fwd*Fb = (Nf*Nb)/(Df*Db)
  // 1+ = (Df*Db + Nf*Nb)/(Df*Db)
  // feedback = (Nf/Df) / ((Df*Db + Nf*Nb)/(Df*Db)) = Nf*Db / (Df*Db + Nf*Nb)
  const num = polyMul(numFwd, denFb);
  const den = polyAdd(polyMul(denFwd, denFb), polyMul(numFwd, numFb));
  return { num, den };
}

function companionEigenvaluesRealPoly(den) {
  // den must be an array like [a0, a1, ..., an] (descending powers)
  if (!Array.isArray(den) || den.length < 2) {
    // Return empty list instead of crashing the whole app
    return [];
  }
  if (!Number.isFinite(den[0]) || Math.abs(den[0]) < 1e-14) {
    return [];
  }

  const a0 = den[0];
  const d = den.map(v => v / a0);
  const n = d.length - 1;
  if (n <= 0) return [];

  const C = numeric.rep([n, n], 0);
  for (let i = 0; i < n - 1; i++) C[i][i + 1] = 1;

  // correct ordering: [-an, ..., -a1]
  for (let j = 0; j < n; j++) C[n - 1][j] = -d[n - j];

  const eig = numeric.eig(C);
  if (!eig || !eig.lambda || !eig.lambda.x || !eig.lambda.y) return [];

  const re = eig.lambda.x;
  const im = eig.lambda.y;

  const out = [];
  for (let i = 0; i < re.length; i++) out.push({ re: re[i], im: im[i] });
  return out;
}


// Random normal (Box–Muller)
let rngSeed = 1234567;

// Axis/time scaling: k in {-2,-1,0,1,2}. Frequency range shifts by ×10^k,
// simulation horizon and timestep shift by ÷10^k.
let axisScaleExp = 0;
const AXIS_SCALE_EXP_MIN = -2;
const AXIS_SCALE_EXP_MAX = 2;

function getAxisScale() {
  return Math.pow(10, axisScaleExp);
}

function setAxisScaleExp(k) {
  axisScaleExp = Math.max(AXIS_SCALE_EXP_MIN, Math.min(AXIS_SCALE_EXP_MAX, k));
  cachedSignals = buildSignals(); // rebuild time signals at new scale
  updateScaleUI();
  renderAll();
}

function fmtPow10Exp(exp) {
  if (exp === 0) return "1";
  return `1e${exp}`;
}

function fmtNice(x) {
  if (!Number.isFinite(x)) return "—";
  // readable, without long tails
  if (Math.abs(x) >= 100) return x.toFixed(0);
  if (Math.abs(x) >= 10) return x.toFixed(1);
  if (Math.abs(x) >= 1) return x.toFixed(2);
  return x.toFixed(3);
}

function updateScaleUI() {
  const scaleEl = document.getElementById("scaleVal");
  const wRangeEl = document.getElementById("wRangeText");
  const tRangeEl = document.getElementById("tRangeText");
  const dtEl = document.getElementById("dtText");
  if (!scaleEl || !wRangeEl || !tRangeEl || !dtEl) return;

  const wMinExp = -3 + axisScaleExp;
  const wMaxExp =  2 + axisScaleExp;
  const scale = getAxisScale();

  scaleEl.textContent = `×${fmtPow10Exp(axisScaleExp)}`;

  wRangeEl.textContent = `${fmtPow10Exp(wMinExp)} … ${fmtPow10Exp(wMaxExp)}`;

  const Tend = 60 / scale;
  const dt = 0.025 / scale;
  tRangeEl.textContent = `0 … ${fmtNice(Tend)}`;
  dtEl.textContent = fmtNice(dt);
}

function randUniform() {
  // xorshift32
  let x = rngSeed | 0;
  x ^= x << 13; x ^= x >>> 17; x ^= x << 5;
  rngSeed = x;
  return ((x >>> 0) / 4294967296);
}
function randn() {
  let u1 = Math.max(randUniform(), 1e-12);
  let u2 = randUniform();
  return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
}

function clamp(x, lo, hi) { return Math.min(hi, Math.max(lo, x)); }
function safeLog10(x) { return Math.log(x) / Math.LN10; }
function setSliderFromValue(slider, value) { slider.value = String(safeLog10(value)); }


// UI hooks
const el = id => document.getElementById(id);
const kSlider = el("kSlider");
const aSlider = el("aSlider");
const bSlider = el("bSlider");
const nSlider = el("nSlider");
const gNumEl = el("gNum");
const gDenEl = el("gDen");
const kInput = el("kInput");
const aInput = el("aInput");
const bInput = el("bInput");
const nInput = el("nInput");

const scaleDownBtn = el("scaleDownBtn");
const scaleUpBtn = el("scaleUpBtn");
const scaleResetBtn = el("scaleResetBtn");
const challengeSelect = el("challengeSelect");
const loadChallengeBtn = el("loadChallengeBtn");


let override = { K: false, a: false, b: false, N: false };

// -----------------------------
// Challenges (Phase 1.5 + Phase 2)
// -----------------------------

// Phase 2 time-domain mask for y(t).
// Use `lo: null` or `hi: null` to indicate an unbounded side.
// User-provided mask spec:
//  - Max overshoot: y < 1.2 for 0 < t < 4
//  - Settling band: 0.95 < y < 1.05 for 5 < t < 20
//  - Settling band: 0.95 < y < 1.05 for 25 < t < 40
const DEFAULT_Y_MASK = [
  { t0: 0,  t1: 5,  lo: null,  hi: 1.2  },
  { t0: 5,  t1: 20, lo: 0.95, hi: 1.05 },
  { t0: 20, t1: 40, lo: 0.7, hi: 1.3 },
];

const CHALLENGES = {
  c1: {
    id: "c1",
    name: "Robust & calm",
    goal: "Achieve robust stability margins with a smooth step response.",
    plantNum: [1],
    plantDen: [1, 3, 2, 0],   // 1/(s(s+1)(s+2))
    defaults: { K: 1, a: 0.01, b: 1, N: 2 },
    scaleExp: 0,
    req: { AmMin: 2.0, PmMin: 45.0 },
    yMask: DEFAULT_Y_MASK,
    maskTol: 1e-3,
    maskMaxViolations: 0,
    hint: "Add phase lead and increase cross-over frequency by (K,b,N). Tune integral action by (a)"
  },
  c2: {
    id: "c2",
    name: "Fast but safe",
    goal: "Achieve higher bandwidth (faster response) while maintaining safe margins.",
    plantNum: [1],
    plantDen: [1, 3, 2, 0],
    defaults: { K: 2, a: 0.02, b: 2, N: 3 },
    scaleExp: 1,
    req: { AmMin: 1.6, PmMin: 35.0 },
    yMask: DEFAULT_Y_MASK,
    maskTol: 1e-3,
    maskMaxViolations: 0,
    hint: "Push bandwidth higher but keep it stable. Increase K, then recover phase margin using lead (b,N). Watch noise sensitivity."
  },
  c3: {
    id: "c3",
    name: "Disturbance fighter",
    goal: "Reject load disturbances quickly without sacrificing robustness.",
    plantNum: [1],
    plantDen: [1, 3, 2, 0],
    defaults: { K: 1.5, a: 0.01, b: 1, N: 2.5 },
    scaleExp: 0,
    req: { AmMin: 2.0, PmMin: 40.0 },
    yMask: DEFAULT_Y_MASK,
    maskTol: 1e-3,
    maskMaxViolations: 0,
    hint: "Focus on rejecting the load step. Keep integral action effective, then add lead to preserve margin."
  }
};

function getSelectedChallenge() {
  if (!challengeSelect) return null;
  const id = challengeSelect.value;
  if (id === "custom") return null;
  return CHALLENGES[id] || null;
}

function setInputsFromControllerDefaults(def) {
  // Sliders take control
  override.K = override.a = override.b = override.N = false;

  // Set sliders (log-scale). Our presets are within slider ranges.
  kSlider.value = String(Math.log10(def.K));
  aSlider.value = String(Math.log10(def.a));
  bSlider.value = String(Math.log10(def.b));
  nSlider.value = String(Math.log10(def.N));

  // Sync numeric boxes for clarity
  kInput.value = def.K.toFixed(6);
  aInput.value = def.a.toFixed(6);
  bInput.value = def.b.toFixed(6);
  nInput.value = def.N.toFixed(6);
}

function applyChallenge(ch) {
  if (!ch) return;

  // Plant
  gNumEl.value = ch.plantNum.join(" ");
  gDenEl.value = ch.plantDen.join(" ");

  // Scale (also rebuild time signals)
  axisScaleExp = ch.scaleExp ?? 0;
  updateScaleUI();
  cachedSignals = buildSignals();

  // Controller defaults
  setInputsFromControllerDefaults(ch.defaults);

  // Hint text
  const hintEl = el("challengeHint");
  if (hintEl) hintEl.textContent = ch.hint;

  renderAll();
}

function statusBadge(ok) { return ok ? "✅" : "❌"; }


function updateChallengeGoalUI(ch) {
  const goalEl = el("goalText");
  const targetEl = el("targetText");
  if (!goalEl || !targetEl) return;

  if (!ch) {
    goalEl.textContent = "—";
    targetEl.textContent = "—";
    return;
  }
  goalEl.textContent = ch.goal || "—";
  const hasMask = Array.isArray(ch.yMask) && ch.yMask.length > 0;
  targetEl.textContent = `Am ≥ ${ch.req.AmMin}, Phim ≥ ${ch.req.PmMin}°` + (hasMask ? `, mask: y(t)` : "");
}

// Simple score (Phase 1.5): rewards exceeding targets, but only if targets are met.
// Score is shown as — until a challenge is selected; then 0 if not passed; otherwise 100+ bonuses.
function computeChallengeScore(ch, Gm, Pm, maskEval) {
  if (!ch) return null;
  // Gain margin can be Infinity (valid, means no phase crossover)
  const gmOkNumber = (Gm === Infinity) || Number.isFinite(Gm);
  if (!gmOkNumber || !Number.isFinite(Pm)) return 0;

  const amOk = Gm >= ch.req.AmMin;
  const pmOk = Pm >= ch.req.PmMin;
  if (!(amOk && pmOk)) return 0;

  // If a time-domain mask is defined, it must pass too
  const hasMask = Array.isArray(ch.yMask) && ch.yMask.length > 0;
  if (hasMask) {
    const pass = maskEval && maskEval.pass;
    if (!pass) return 0;
  }

  // Bonuses for exceeding margins (kept small and easy to understand)
  // Cap bonus if Gm is infinite
  const gmForBonus = (Gm === Infinity) ? (ch.req.AmMin + 2) : Gm;
  const amBonus = 20 * Math.max(0, Math.min(2, (gmForBonus - ch.req.AmMin)));  // up to +40
  const pmBonus = 1.0 * Math.max(0, Math.min(30, (Pm - ch.req.PmMin))); // up to +30

  return Math.round(100 + amBonus + pmBonus);
}

function updateChallengeStatusUI(ch, Gm, Pm, maskEval) {
  const scoreEl = el("statusScore");
  updateChallengeGoalUI(ch);

  const amEl = el("statusAm");
  const pmEl = el("statusPm");
  const maskEl = el("statusMask");
  if (!amEl || !pmEl) return;

  if (!ch) {
    amEl.textContent = `Am: ${Number.isFinite(Gm) ? Gm.toFixed(2) : (Gm === Infinity ? "∞" : "—")}`;
    pmEl.textContent = `Phim: ${Number.isFinite(Pm) ? Pm.toFixed(1) : "—"}°`;
    if (maskEl) maskEl.textContent = "Mask: —";
    if (scoreEl) scoreEl.textContent = "Score: —";
    return;
  }

  const amOk = ((Gm === Infinity) || Number.isFinite(Gm)) && (Gm >= ch.req.AmMin);
  const pmOk = Number.isFinite(Pm) && Pm >= ch.req.PmMin;

  amEl.textContent = `Am: ${Number.isFinite(Gm) ? Gm.toFixed(2) : (Gm === Infinity ? "∞" : "—")} ${statusBadge(amOk)} (≥ ${ch.req.AmMin})`;
  pmEl.textContent = `Phim: ${Number.isFinite(Pm) ? Pm.toFixed(1) : "—"}° ${statusBadge(pmOk)} (≥ ${ch.req.PmMin}°)`;

  const hasMask = Array.isArray(ch.yMask) && ch.yMask.length > 0;
  if (maskEl) {
    if (!hasMask || !maskEval) {
      maskEl.textContent = "Mask: —";
    } else {
      const ok = !!maskEval.pass;
      maskEl.textContent = `Mask: ${statusBadge(ok)} (${maskEval.violationTime.toFixed(2)} s outside)`;
    }
  }

  const score = computeChallengeScore(ch, Gm, Pm, maskEval);
  if (scoreEl) scoreEl.textContent = (score === null) ? "Score: —" : `Score: ${score}`;
}

// -----------------------------
// Phase 2: time-domain mask
// -----------------------------

function buildYMaskShapes(maskSegments, yMin, yMax) {
  if (!maskSegments || !maskSegments.length) return [];
  const y0Fallback = Number.isFinite(yMin) ? yMin : -1;
  const y1Fallback = Number.isFinite(yMax) ? yMax : 2;
  return maskSegments.map(seg => {
    const y0 = (seg.lo === null || !Number.isFinite(seg.lo)) ? y0Fallback : seg.lo;
    const y1 = (seg.hi === null || !Number.isFinite(seg.hi)) ? y1Fallback : seg.hi;
    return {
      type: "rect",
      xref: "x",
      yref: "y",
      x0: seg.t0,
      x1: seg.t1,
      y0,
      y1,
      line: { width: 0 },
      fillcolor: "rgba(0, 180, 0, 0.12)",
      layer: "below"
    };
  });
}

function maskBoundsAtTime(maskSegments, t) {
  if (!maskSegments) return null;
  for (let i = 0; i < maskSegments.length; i++) {
    const s = maskSegments[i];
    if (t >= s.t0 && t < s.t1) return { lo: s.lo, hi: s.hi };
  }
  return null;
}

function evaluateMask(tArr, yArr, maskSegments, tol = 1e-3, maxViolations = 0) {
  if (!maskSegments || maskSegments.length === 0) {
    return { pass: true, violationTime: 0, violations: 0 };
  }
  let violations = 0;
  let violationTime = 0;

  for (let k = 0; k < tArr.length; k++) {
    const t = tArr[k];
    const b = maskBoundsAtTime(maskSegments, t);
    if (!b) continue;

    const y = yArr[k];
    const loOk = (b.lo === null || !Number.isFinite(b.lo)) ? true : (y >= b.lo - tol);
    const hiOk = (b.hi === null || !Number.isFinite(b.hi)) ? true : (y <= b.hi + tol);
    const bad = !(loOk && hiOk);
    if (bad) {
      violations++;
      const dt = (k === 0) ? 0 : (tArr[k] - tArr[k - 1]);
      violationTime += dt;
    }
  }

  const pass = (violations <= (maxViolations ?? 0));
  return { pass, violationTime, violations };
}



function sliderVal(sl) { return parseFloat(sl.value); }
function pow10(x) { return Math.pow(10, x); }

function getParams() {
  // default: slider values
  let K = pow10(sliderVal(kSlider));
  let a = pow10(sliderVal(aSlider));
  let b = pow10(sliderVal(bSlider));
  let N = pow10(sliderVal(nSlider));

  // only override if that input is "active override"
  if (override.K) {
    const v = Number(kInput.value);
    if (Number.isFinite(v) && v > 0) K = v;
  }
  if (override.a) {
    const v = Number(aInput.value);
    if (Number.isFinite(v) && v > 0) a = v;
  }
  if (override.b) {
    const v = Number(bInput.value);
    if (Number.isFinite(v) && v > 0) b = v;
  }
  if (override.N) {
    const v = Number(nInput.value);
    if (Number.isFinite(v) && v > 0) N = v;
  }

  return { K, a, b, N };
}


function updateParamLabels() {
  const { K, a, b, N } = getParams();

  // Display values (6 decimals)
  el("kVal").textContent = K.toFixed(6);
  el("aVal").textContent = a.toFixed(6);
  el("bVal").textContent = b.toFixed(6);
  el("nVal").textContent = N.toFixed(6);

  // Keep number inputs in sync (avoid overwriting while user is typing)
  if (document.activeElement !== kInput) kInput.value = K.toFixed(6);
  if (document.activeElement !== aInput) aInput.value = a.toFixed(6);
  if (document.activeElement !== bInput) bInput.value = b.toFixed(6);
  if (document.activeElement !== nInput) nInput.value = N.toFixed(6);
}



function makeControllerTF(K, a, b, N) {
  // C(s) = K*(s+a)/s * N*(s+b)/(s+bN)
  // (s+a)/s => num=[1 a], den=[1 0]
  // N*(s+b)/(s+bN) => num=[N N*b], den=[1 b*N]
  const num1 = [1, a];
  const den1 = [1, 0];
  const num2 = [N, N * b];
  const den2 = [1, b * N];

  const num = polyScale(polyMul(num1, num2), K);
  const den = polyMul(den1, den2);
  return { num, den };
}

function buildSignals() {
  const scale = getAxisScale();
  const dt = 0.025 / scale;
  const Tend = 60 / scale;

  const t = [];
  const nSteps = Math.round(Tend / dt);
  for (let i = 0; i <= nSteps; i++) t.push(i * dt);

  const tNoise = 40 / scale;
  const tLoad = 20 / scale;

  const r = t.map(() => 1);
  const n = t.map(tt => (tt > tNoise ? 0.2 * randn() : 0));
  const l = t.map(tt => (tt > tLoad ? -0.5 : 0));
  return { t, r, n, l };
}

let cachedSignals = buildSignals();
updateScaleUI();

function renderAll() {
  updateParamLabels();

  let numG, denG;
  try {
    numG = parseCoeffs(gNumEl.value);
    denG = parseCoeffs(gDenEl.value);
  } catch (e) {
    console.warn("Input error:", e.message);
    const hintEl = el("challengeHint");
    if (hintEl) hintEl.textContent = `Input error: ${e.message}`;
    return;
  }

  const { K, a, b, N } = getParams();
  const Ctf = makeControllerTF(K, a, b, N);

  // Frequency response
  const w = logspace(-3 + axisScaleExp, 2 + axisScaleExp, 200);

  // Open loop L = G*C
  const numL = polyMul(numG, Ctf.num);
  const denL = polyMul(denG, Ctf.den);

  const Ljw = w.map(ww => tfFreqResp(numL, denL, ww));
  const mag = Ljw.map(z => cAbs(z));
  const phase = unwrapPhaseDeg(Ljw.map(z => cArgDeg(z)));

  const { Gm, Pm, Wcg, Wcp } = computeMargins(w, mag, phase);

  const ch = getSelectedChallenge();
  // Update margins immediately (mask evaluated after time simulation)
  updateChallengeStatusUI(ch, Gm, Pm, null);

  // Bode plot (two traces: before=G, after=G*C)
  const Gjw = w.map(ww => tfFreqResp(numG, denG, ww));
  const magG = Gjw.map(z => cAbs(z));
  const phaseG = unwrapPhaseDeg(Gjw.map(z => cArgDeg(z)));

    const bodeData = [
    { x: w, y: magG, name: "G (before)",   mode: "lines", xaxis: "x",  yaxis: "y",
        line: { color: "#1f77b4" } }, // blue

    { x: w, y: mag,  name: "G·C (after)",  mode: "lines", xaxis: "x",  yaxis: "y",
        line: { color: "#ff7f0e" } }, // orange

    { x: w, y: phaseG, name: "Phase G",    mode: "lines", xaxis: "x2", yaxis: "y2", showlegend: false,
        line: { color: "#1f77b4" } }, // blue

    { x: w, y: phase,  name: "Phase G·C",  mode: "lines", xaxis: "x2", yaxis: "y2", showlegend: false,
        line: { color: "#ff7f0e" } }, // orange
    ];


  const bodeLayout = {
    margin: { l: 55, r: 10, t: 60, b: 45 },
    title: {
      text: `Bode diagram<br><span style="font-size:11px;">Gain margin Am = ${Number.isFinite(Gm) ? Gm.toFixed(2) : "—"} at ω = ${Number.isFinite(Wcg) ? Wcg.toFixed(2) : "—"} &nbsp; | &nbsp; Phase margin Phim = ${Number.isFinite(Pm) ? Pm.toFixed(2) : "—"} at ω = ${Number.isFinite(Wcp) ? Wcp.toFixed(2) : "—"}</span>`,
      font: { size: 14 }                      // smaller title
    },
    legend: { orientation: "h" },
    grid: { rows: 2, columns: 1, pattern: "independent" },
    xaxis: { type: "log", title: "ω (rad/s)" },
    yaxis: { type: "log", title: "|L(jω)|", range: [-3, 3] },
    xaxis2: { type: "log", title: "ω (rad/s)" },
    yaxis2: { title: "∠L(jω) (deg)", range: [-300, 0] },
    shapes: [
      // |L|=1 line
      { type: "line", xref: "x", yref: "y", x0: w[0], x1: w[w.length - 1], y0: 1, y1: 1, line: { width: 2 } },
      // phase=-180 line
      { type: "line", xref: "x2", yref: "y2", x0: w[0], x1: w[w.length - 1], y0: -180, y1: -180, line: { width: 2 } },
    ]
  };

  Plotly.react("bodePlot", bodeData, bodeLayout, { responsive: true });

  // Nyquist: plot L(jw) and mirrored branch like MATLAB nyquist does
  const xNy = Ljw.map(z => z.re);
  const yNy = Ljw.map(z => z.im);
  const xNyMir = xNy.concat(xNy.slice().reverse());
  const yNyMir = yNy.concat(yNy.slice().reverse().map(v => -v));

  const xNy0 = Gjw.map(z => z.re);
  const yNy0 = Gjw.map(z => z.im);
  const xNy0Mir = xNy0.concat(xNy0.slice().reverse());
  const yNy0Mir = yNy0.concat(yNy0.slice().reverse().map(v => -v));

  const nyData = [
    { x: xNy0Mir, y: yNy0Mir, mode: "lines", name: "G (before)" },
    { x: xNyMir,  y: yNyMir,  mode: "lines", name: "G·C (after)" },
    { x: [-1], y: [0], mode: "markers", name: "-1", marker: { size: 10 } }
  ];
  const nyLayout = {
    margin: { l: 55, r: 10, t: 40, b: 45 },
    title: "Nyquist diagram",
    xaxis: { title: "Re", range: [-3, 3] },
    yaxis: { title: "Im", range: [-3, 3], scaleanchor: "x", scaleratio: 1 },
    legend: { orientation: "h" }
  };
  Plotly.react("nyquistPlot", nyData, nyLayout, { responsive: true });

  // Time domain
  const { t, r, n: noise, l } = cachedSignals;
  const rMinusN = r.map((v, i) => v - noise[i]);

  // Pre-comp (Gc = 1)
  const oneNum = [1], oneDen = [1];

  const T1 = feedbackTF(polyMul(numG, oneNum), polyMul(denG, oneDen), [1], [1]);     // feedback(G*1,1) = G/(1+G)
  const G_over_1pG = feedbackTF(numG, denG, [1], [1]);                               // feedback(G,1)
  const u_from_r_1 = feedbackTF(oneNum, oneDen, numG, denG);                         // feedback(1,G)

  const y1a = simulateTF(T1.num, T1.den, rMinusN, t);
  const y1b = simulateTF(feedbackTF(numG, denG, [1], [1]).num, feedbackTF(numG, denG, [1], [1]).den, l, t); // actually lsim(feedback(G,1), l)
  const y1 = y1a.map((v, i) => v + y1b[i]);

  const u1a = simulateTF(u_from_r_1.num, u_from_r_1.den, rMinusN, t);
  const u1b = simulateTF(T1.num, T1.den, l, t);
  const u1 = u1a.map((v, i) => v + u1b[i]);

  // Compensated
  // T = feedback(G*C,1)
  const GC_num = polyMul(numG, Ctf.num);
  const GC_den = polyMul(denG, Ctf.den);
  const T = feedbackTF(GC_num, GC_den, [1], [1]);

  // G/(1+GC) = feedback(G, C)
  const Gdist = feedbackTF(numG, denG, Ctf.num, Ctf.den);

  // C/(1+GC) = feedback(C, G)
  const Cu = feedbackTF(Ctf.num, Ctf.den, numG, denG);

  const ya = simulateTF(T.num, T.den, rMinusN, t);
  const yb = simulateTF(Gdist.num, Gdist.den, l, t);
  const y = ya.map((v, i) => v + yb[i]);

  const ua = simulateTF(Cu.num, Cu.den, rMinusN, t);
  const ub = simulateTF(T.num, T.den, l, t);
  const u = ua.map((v, i) => v + ub[i]);

  // Phase 2: evaluate y(t) mask (after response) and update challenge status/score
  const hasMask = (ch && Array.isArray(ch.yMask) && ch.yMask.length > 0);
  const maskEval = hasMask ? evaluateMask(t, y, ch.yMask, ch.maskTol ?? 1e-3, ch.maskMaxViolations ?? 0) : null;
  updateChallengeStatusUI(ch, Gm, Pm, maskEval);

  // y-axis autoscale: include response traces and finite mask bounds (ignore null/unbounded)
  const yVals = y.concat(y1).filter(Number.isFinite);
  let yMin = yVals.length ? Math.min(...yVals) : 0;
  let yMax = yVals.length ? Math.max(...yVals) : 1.5;

  if (hasMask) {
    const loVals = ch.yMask.map(s => s.lo).filter(Number.isFinite);
    const hiVals = ch.yMask.map(s => s.hi).filter(Number.isFinite);
    if (loVals.length) yMin = Math.min(yMin, Math.min(...loVals));
    if (hiVals.length) yMax = Math.max(yMax, Math.max(...hiVals));
  }

  // Add a little padding
  const pad = 0.08 * Math.max(1e-6, (yMax - yMin));
  yMin -= pad;
  yMax += pad;

  const yMaskShapes = hasMask ? buildYMaskShapes(ch.yMask, yMin, yMax) : [];

  Plotly.react("yPlot", [
    // draw after first (orange)
    { x: t, y: y,  mode: "lines", name: "y after",  line: { color: "#ff7f0e" } },

    // draw before last (blue) so it sits on top
    { x: t, y: y1, mode: "lines", name: "y before", line: { color: "#1f77b4" } },

    // reference line
    { x: [t[0], t[t.length - 1]], y: [1, 1],
      mode: "lines",
      name: "reference",
      line: { dash: "dash", color: "#444" } },
  ], {
    margin: { l: 55, r: 10, t: 20, b: 45 },
    yaxis: { title: "output y", range: [yMin, yMax] },
    xaxis: { title: "time (s)", range: [0, t[t.length - 1]] },
    legend: { orientation: "h" },
    shapes: yMaskShapes
  }, { responsive: true });

  Plotly.react("uPlot", [
    { x: t, y: u,  mode: "lines", name: "u after",  line: { color: "#ff7f0e" } },
    { x: t, y: u1, mode: "lines", name: "u before", line: { color: "#1f77b4" } },
  ], {
    margin: { l: 55, r: 10, t: 20, b: 45 },
    yaxis: { title: "control u" },
    xaxis: { title: "time (s)", range: [0, t[t.length - 1]] },
    legend: { orientation: "h" }
  }, { responsive: true });

  // Poles (before = feedback(G,1), after = feedback(G*C,1))
 // Poles (before = feedback(G,1), after = feedback(G*C,1))
let polesBefore = [];
let polesAfter = [];
let poleDebug = "";

try {
  // Closed-loop denominators via the same helper used elsewhere
  const beforeCL = feedbackTF(numG, denG, [1], [1]);   // G/(1+G)
  const afterCL  = feedbackTF(numL, denL, [1], [1]);   // GC/(1+GC)

  //polesBefore = companionEigenvaluesRealPoly(beforeCL.den);
  //polesAfter  = companionEigenvaluesRealPoly(afterCL.den);

  polesBefore = durandKernerRootsRealPoly(beforeCL.den);
  polesAfter  = durandKernerRootsRealPoly(afterCL.den);

  poleDebug = `
    beforeCL.den length: ${beforeCL.den?.length ?? "undef"}<br>
    afterCL.den length: ${afterCL.den?.length ?? "undef"}<br>
    polesBefore: ${polesBefore.length}, polesAfter: ${polesAfter.length}<br>
  `;
} catch (e) {
  poleDebug = `<span style="color:#b00;"><b>Pole calc failed:</b> ${e.message}</span>`;
}

// Auto axis limits so poles are always visible
const allPoles = polesBefore.concat(polesAfter);
const reVals = allPoles.map(p => p.re).filter(Number.isFinite);
const imVals = allPoles.map(p => p.im).filter(Number.isFinite);

let reMin = -3, reMax = 1, imMax = 3;
if (reVals.length) {
  const rmin = Math.min(...reVals);
  const rmax = Math.max(...reVals);
  const pad = 0.2 * Math.max(1, rmax - rmin);
  reMin = rmin - pad;
  reMax = rmax + pad;
}
if (imVals.length) {
  const im = Math.max(...imVals.map(v => Math.abs(v)));
  imMax = Math.max(1, im + 0.5);
}

Plotly.react("polePlot", [
  {
    x: polesBefore.map(p => p.re),
    y: polesBefore.map(p => p.im),
    mode: "markers",
    name: "before",
    marker: { symbol: "x", size: 10 }
  },
  {
    x: polesAfter.map(p => p.re),
    y: polesAfter.map(p => p.im),
    mode: "markers",
    name: "after",
    marker: { symbol: "x", size: 10 }
  },
  // axes lines
  { x: [0, 0], y: [-imMax, imMax], mode: "lines", name: "Re=0", line: { dash: "dash" }, showlegend: false },
  { x: [reMin, reMax], y: [0, 0], mode: "lines", name: "Im=0", line: { dash: "dash" }, showlegend: false },
], {
  margin: { l: 55, r: 10, t: 40, b: 45 },
  title: "Poles (blue: before, red: after)",
  xaxis: { title: "Re", range: [reMin, reMax] },
  yaxis: { title: "Im", range: [-imMax, imMax], scaleanchor: "x", scaleratio: 1 },
  legend: { orientation: "h" }
}, { responsive: true });

// Optional: show poleDebug in infoBox (or console)
// console.log(poleDebug);


  // Info box
 // el("infoBox").innerHTML = `
 //   <div style="font-weight:800; font-size:16px; margin-bottom:8px;">Current setup</div>
 //   <div><b>G(s)</b> = [${numG.join(", ")}] / [${denG.join(", ")}]</div>
 //   <div style="margin-top:8px;"><b>K</b> = ${K.toPrecision(4)}, <b>a</b> = ${a.toPrecision(4)}, <b>b</b> = ${b.toPrecision(4)}, <b>N</b> = ${N.toPrecision(4)}</div>
 //   <div style="margin-top:8px;"><b>Margins (approx)</b></div>
 //   <div>Am = ${Number.isFinite(Gm) ? Gm.toFixed(3) : "—"} at ωcg = ${Number.isFinite(Wcg) ? Wcg.toFixed(3) : "—"}</div>
 //   <div>Phim = ${Number.isFinite(Pm) ? Pm.toFixed(2) : "—"}° at ωcp = ${Number.isFinite(Wcp) ? Wcp.toFixed(3) : "—"}</div>
 //   <div style="margin-top:10px;" class="small">
 //     Note: margins are computed on the sampled ω-grid with log-interpolation.
 //   </div>
 // `;
}

// Events
kSlider.addEventListener("input", () => { override.K = false; renderAll(); });
aSlider.addEventListener("input", () => { override.a = false; renderAll(); });
bSlider.addEventListener("input", () => { override.b = false; renderAll(); });
nSlider.addEventListener("input", () => { override.N = false; renderAll(); });

[gNumEl, gDenEl].forEach(inp => inp.addEventListener("change", renderAll));


if (scaleDownBtn) scaleDownBtn.addEventListener("click", () => setAxisScaleExp(axisScaleExp - 1));
if (scaleUpBtn) scaleUpBtn.addEventListener("click", () => setAxisScaleExp(axisScaleExp + 1));
if (scaleResetBtn) scaleResetBtn.addEventListener("click", () => setAxisScaleExp(0));


if (loadChallengeBtn) loadChallengeBtn.addEventListener("click", () => {
  const ch = getSelectedChallenge();
  if (ch) applyChallenge(ch);
});

if (challengeSelect) challengeSelect.addEventListener("change", () => {
  const ch = getSelectedChallenge();
  const hintEl = el("challengeHint");
  if (hintEl) hintEl.textContent = ch ? ch.hint : "Custom mode: use your own plant and goals.";
  updateChallengeGoalUI(ch);
  renderAll();
});

el("reseedBtn").addEventListener("click", () => {
  rngSeed = (Date.now() & 0xffffffff) ^ 0x9e3779b9;
  cachedSignals = buildSignals();
  renderAll();
});

el("resetBtn").addEventListener("click", () => {
  axisScaleExp = 0;
  updateScaleUI();
  gNumEl.value = "1";
  gDenEl.value = "1 3 2 0";
  kSlider.value = "0";        // K=1
  aSlider.value = "-2";       // a=0.01
  bSlider.value = "0";        // b=1
  nSlider.value = "0.30103";  // N≈2
  rngSeed = 1234567;
  cachedSignals = buildSignals();
  renderAll();
});

function hookNumericInput(inputEl, sliderEl, which, sliderMinLog, sliderMaxLog) {
  function apply() {
    let v = Number(inputEl.value);
    if (!Number.isFinite(v) || v <= 0) return;

    // typing activates override for this parameter
    override[which] = true;

    // move slider only if in range (no clamp)
    const logv = Math.log10(v);
    if (logv >= sliderMinLog && logv <= sliderMaxLog) {
      sliderEl.value = String(logv);
    }

    renderAll();
  }

  inputEl.addEventListener("change", apply);
  inputEl.addEventListener("keydown", (e) => {
    if (e.key === "Enter") apply();
  });
}


// Use the sliders' log ranges directly:
hookNumericInput(kInput, kSlider, "K", -1, 1);
hookNumericInput(aInput, aSlider, "a", -2, 0);
hookNumericInput(bInput, bSlider, "b", -1, 1);
hookNumericInput(nInput, nSlider, "N", 0, 1.30103);

// Initial render
updateScaleUI();
updateParamLabels();
renderAll();
