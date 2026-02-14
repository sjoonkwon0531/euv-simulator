import { useState, useCallback, useMemo, useRef, useEffect } from "react";
import * as math from "mathjs";
import _ from "lodash";
import {
  LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend,
  ResponsiveContainer, ScatterChart, Scatter, AreaChart, Area,
  RadarChart, PolarGrid, PolarAngleAxis, PolarRadiusAxis, Radar,
  ComposedChart, ReferenceLine, BarChart, Bar
} from "recharts";

// ════════════════════════════════════════════════════════════════
// CONSTANTS & PHYSICS ENGINE
// ════════════════════════════════════════════════════════════════

const C = {
  h: 6.626e-34, c: 2.998e8, eV: 1.602e-19, kB: 1.381e-23,
  EUV_wl: 13.5e-9, EUV_eV: 91.8,
  Sn_sigma: 2.8e-20, // cm² per Sn atom at 13.5nm
};
const EUV_E_J = C.h * C.c / C.EUV_wl;

function photonDensity(dose) { return dose * 10 / EUV_E_J; }
function absCoeff(nSn, vol) { return nSn * C.Sn_sigma * 1e14 / vol; }
function photonsPerMol(dose, mu, sz) {
  const pd = photonDensity(dose) * 1e-18;
  return pd * sz * sz * (1 - Math.exp(-mu * sz));
}
function binomPMF(k, n, p) {
  if (n < 0 || k < 0 || k > n) return 0;
  if (p <= 0) return k === 0 ? 1 : 0;
  if (p >= 1) return k === n ? 1 : 0;
  return math.combinations(n, k) * Math.pow(p, k) * Math.pow(1 - p, n - k);
}
function cumBinom(thresh, n, p) {
  let s = 0;
  for (let k = thresh; k <= n; k++) s += binomPMF(k, n, p);
  return Math.min(1, Math.max(0, s));
}
function agingModel(t, T, RH, params) {
  const Tk = T + 273.15, kBT = C.kB * Tk / C.eV;
  const rH = (params.k0h || 1e-3) * Math.exp(-(params.Eah || 0.65) / kBT) * (RH / 50);
  const rA = (params.k0a || 5e-4) * Math.exp(-(params.Eaa || 0.45) / kBT);
  const fH = Math.exp(-rH * t), fA = Math.exp(-rA * t);
  return {
    intact: fH * fA,
    sizeFactor: 1 + 0.15 * (1 - fA),
    absFactor: fH * 0.95 + 0.05,
    odIncrease: (1 - fH) * 0.3 + (1 - fA) * 0.2,
  };
}

function runModel(inp) {
  const { nSn, molSz, ringD, filmTh, dose, NA, hp, seB, chemB, thRE, thSCF, agT, temp, hum, TON } = inp;
  const wl = 13.5, vol = (4/3) * Math.PI * Math.pow(ringD / 2, 3);
  const mu = absCoeff(nSn, vol);
  const ppm = photonsPerMol(dose, mu, molSz);
  const blur = Math.sqrt(seB ** 2 + chemB ** 2);
  const nSites = Math.max(1, Math.round(vol * 8));
  const rpm = ppm * TON;
  const rProb = Math.min(0.99, rpm / nSites);
  const pMSS_e = cumBinom(thRE, nSites, rProb);
  const pMSS_u = cumBinom(thRE, nSites, rProb * 0.05);
  const pSC_e = cumBinom(thSCF, 27, pMSS_e);
  const pSC_u = cumBinom(thSCF, 27, pMSS_u);
  const nMolTh = Math.round(filmTh / molSz);
  const k1 = hp * NA / wl;
  const imgC = Math.min(1, Math.max(0, 1 - 0.5 * Math.exp(-3 * (k1 - 0.25))));

  // Spatial profile
  const nP = 50, drProf = [];
  for (let i = 0; i < nP; i++) {
    const x = (i / (nP - 1) - 0.5) * 2 * hp;
    const img = 0.5 * (1 + imgC * Math.cos(Math.PI * x / hp));
    const lrp = Math.min(0.99, rProb * img / 0.5);
    const lMSS = cumBinom(thRE, nSites, lrp);
    const lSC = cumBinom(thSCF, 27, lMSS);
    const dr = nMolTh * lSC, drV = nMolTh * lSC * (1 - lSC);
    drProf.push({ x: +x.toFixed(2), dr: +dr.toFixed(3), drV: +drV.toFixed(3), img: +img.toFixed(3) });
  }

  const eIdx = Math.floor(nP / 4), dx = drProf[1].x - drProf[0].x;
  const drSlope = Math.abs((drProf[eIdx + 1].dr - drProf[eIdx - 1].dr) / (2 * dx));
  const drVar = drProf[eIdx].drV;
  const drE = drProf[eIdx].dr;
  const uLS = drE > 0 ? (drSlope / drE) * (wl / NA) : 0;

  const calcLER = (sl, v) => {
    if (sl <= 0) return hp * 0.3;
    const base = Math.sqrt(v) / sl;
    const corr = k1 < 0.35 ? 1 + 2 * (0.35 - k1) : 1;
    return Math.min(base * corr * 3, hp * 0.3);
  };

  const ler = calcLER(drSlope, drVar), lwr = ler * Math.SQRT2;
  const defP = Math.max(1 - pSC_e, pSC_u);

  // Aging
  const ag = agingModel(agT, temp, hum, {});
  const agMu = mu * ag.absFactor, agSz = molSz * ag.sizeFactor;
  const agPPM = photonsPerMol(dose, agMu, agSz);
  const agRP = Math.min(0.99, agPPM * TON / nSites);
  const agMSS = cumBinom(thRE, nSites, agRP);
  const agSC = cumBinom(thSCF, 27, agMSS);
  const agOD = 0.05 + ag.odIncrease;
  const agDRV = drVar * (1 + agOD * 5);
  const agLER = calcLER(drSlope * ag.absFactor, agDRV);
  const doseCor = ag.absFactor > 0.01 ? 1 / ag.absFactor : 100;

  // Aging timeline
  const tPts = [0, 1, 2, 4, 8, 12, 24, 48, 72, 96, 120, 168, 336, 720];
  const agTL = tPts.map(t => {
    const a = agingModel(t, temp, hum, {});
    const aM = mu * a.absFactor, aS = molSz * a.sizeFactor;
    const ap = photonsPerMol(dose, aM, aS);
    const ar = Math.min(0.99, ap * TON / nSites);
    const am = cumBinom(thRE, nSites, ar);
    const as2 = cumBinom(thSCF, 27, am);
    const aod = 0.05 + a.odIncrease;
    const adv = drVar * (1 + aod * 5);
    const al = calcLER(drSlope * a.absFactor, adv);
    return {
      time: t, tLabel: t < 24 ? `${t}h` : `${Math.round(t / 24)}d`,
      intact: +(a.intact * 100).toFixed(1), ler: +al.toFixed(3), lwr: +(al * Math.SQRT2).toFixed(3),
      absPct: +(a.absFactor * 100).toFixed(1), doseCor: +(1 / Math.max(0.01, a.absFactor)).toFixed(2),
      szInc: +(a.sizeFactor * 100 - 100).toFixed(1),
    };
  });

  // Dose sweep
  const dSweep = [];
  for (let d = 5; d <= 120; d += 2) {
    const dp = photonsPerMol(d, mu, molSz);
    const rp2 = Math.min(0.99, dp * TON / nSites);
    const m2 = cumBinom(thRE, nSites, rp2);
    const s2 = cumBinom(thSCF, 27, m2);
    const dr2 = nMolTh * s2, dv2 = nMolTh * s2 * (1 - s2);
    const sl2 = drSlope * (d / dose);
    const l2 = calcLER(sl2, dv2);
    const c2 = s2 - cumBinom(thSCF, 27, cumBinom(thRE, nSites, rp2 * 0.05));
    dSweep.push({ dose: d, ler: +l2.toFixed(3), lwr: +(l2 * Math.SQRT2).toFixed(3), contrast: +c2.toFixed(3) });
  }

  // Radar
  const radar = [
    { m: "Resolution", v: +(Math.min(1, k1 / 0.5) * 100).toFixed(0) },
    { m: "LER margin", v: +(Math.max(0, 1 - ler / 3) * 100).toFixed(0) },
    { m: "LWR margin", v: +(Math.max(0, 1 - lwr / 4.5) * 100).toFixed(0) },
    { m: "Defect", v: +(Math.max(0, -Math.log10(Math.max(1e-12, defP)) / 12) * 100).toFixed(0) },
    { m: "Sensitivity", v: +(Math.max(0, 1 - dose / 100) * 100).toFixed(0) },
    { m: "Stability", v: +(ag.intact * 100).toFixed(0) },
  ];

  return { mu, ppm, rpm, rProb, pMSS_e, pMSS_u, pSC_e, pSC_u, uLS, k1, imgC, ler, lwr, defP, nMolTh, drSlope, drVar,
    ag, agLER, agLWR: agLER * Math.SQRT2, doseCor, drProf, agTL, dSweep, radar, blur, nSites };
}

// ════════════════════════════════════════════════════════════════
// IMAGE & DATA ANALYSIS ENGINE
// ════════════════════════════════════════════════════════════════

function analyzeImageForLER(imageData, width, height) {
  // Edge detection via Sobel-like gradient on grayscale
  const gray = new Float32Array(width * height);
  for (let i = 0; i < width * height; i++) {
    gray[i] = (imageData[i * 4] * 0.299 + imageData[i * 4 + 1] * 0.587 + imageData[i * 4 + 2] * 0.114) / 255;
  }

  // Otsu threshold
  let hist = new Array(256).fill(0);
  for (let i = 0; i < gray.length; i++) hist[Math.round(gray[i] * 255)]++;
  const total = gray.length;
  let sumAll = 0;
  for (let i = 0; i < 256; i++) sumAll += i * hist[i];
  let sumB = 0, wB = 0, maxVar = 0, threshold = 128;
  for (let i = 0; i < 256; i++) {
    wB += hist[i];
    if (wB === 0) continue;
    const wF = total - wB;
    if (wF === 0) break;
    sumB += i * hist[i];
    const mB = sumB / wB, mF = (sumAll - sumB) / wF;
    const v = wB * wF * (mB - mF) * (mB - mF);
    if (v > maxVar) { maxVar = v; threshold = i; }
  }
  const th = threshold / 255;

  // Find edges per row
  const edgePositions = [];
  for (let y = 0; y < height; y++) {
    const edges = [];
    for (let x = 1; x < width; x++) {
      const prev = gray[y * width + x - 1], curr = gray[y * width + x];
      if ((prev < th && curr >= th) || (prev >= th && curr < th)) {
        edges.push(x);
      }
    }
    if (edges.length >= 2) edgePositions.push(edges);
  }

  if (edgePositions.length < 10) {
    return { success: false, message: "Could not detect sufficient edges. Ensure image shows clear line/space pattern." };
  }

  // Calculate edge roughness statistics
  const leftEdges = edgePositions.map(e => e[0]);
  const rightEdges = edgePositions.map(e => e[e.length > 1 ? 1 : 0]);

  const stats = (arr) => {
    const mean = arr.reduce((a, b) => a + b, 0) / arr.length;
    const std = Math.sqrt(arr.reduce((s, v) => s + (v - mean) ** 2, 0) / arr.length);
    return { mean, std, sigma3: std * 3 };
  };

  const leftStats = stats(leftEdges);
  const rightStats = stats(rightEdges);

  // Width statistics
  const widths = edgePositions.filter(e => e.length >= 2).map(e => e[1] - e[0]);
  const widthStats = stats(widths);

  // Power spectrum (simplified)
  const N = Math.min(256, leftEdges.length);
  const edgeSample = leftEdges.slice(0, N);
  const mean = edgeSample.reduce((a, b) => a + b, 0) / N;
  const centered = edgeSample.map(v => v - mean);
  const psd = [];
  for (let f = 1; f <= N / 2; f++) {
    let re = 0, im = 0;
    for (let n = 0; n < N; n++) {
      re += centered[n] * Math.cos(2 * Math.PI * f * n / N);
      im -= centered[n] * Math.sin(2 * Math.PI * f * n / N);
    }
    psd.push({ freq: f / N, power: (re * re + im * im) / (N * N) });
  }

  // Edge profile data for charting
  const edgeProfile = leftEdges.slice(0, 200).map((v, i) => ({ y: i, x: v, xR: rightEdges[i] || v + widthStats.mean }));

  return {
    success: true,
    lerLeft: leftStats.sigma3,
    lerRight: rightStats.sigma3,
    lerAvg: (leftStats.sigma3 + rightStats.sigma3) / 2,
    lwr: widthStats.sigma3,
    cdMean: widthStats.mean,
    cdStd: widthStats.std,
    threshold: th,
    nEdges: edgePositions.length,
    psd: psd.slice(0, 40),
    edgeProfile,
    imageWidth: width,
    imageHeight: height,
  };
}

function parseCSV(text) {
  const lines = text.trim().split("\n").filter(l => l.trim());
  if (lines.length < 2) return null;
  const headers = lines[0].split(/[,\t;]/).map(h => h.trim().toLowerCase());
  const data = [];
  for (let i = 1; i < lines.length; i++) {
    const vals = lines[i].split(/[,\t;]/).map(v => v.trim());
    const row = {};
    headers.forEach((h, j) => {
      row[h] = isNaN(parseFloat(vals[j])) ? vals[j] : parseFloat(vals[j]);
    });
    data.push(row);
  }
  return { headers, data };
}

function fitDoseLER(data) {
  // Fit LER = A / sqrt(dose) + B  (photon shot noise model)
  // Or more generally LER = A * dose^(-alpha) + B
  const doseKey = data.headers.find(h => h.includes("dose") || h.includes("energy"));
  const lerKey = data.headers.find(h => h.includes("ler") || h.includes("roughness") || h.includes("lwr") || h.includes("sigma"));
  if (!doseKey || !lerKey) return { success: false, message: `Headers need 'dose' and 'ler/roughness' columns. Found: ${data.headers.join(", ")}` };

  const pts = data.data.filter(r => r[doseKey] > 0 && r[lerKey] > 0).map(r => ({ dose: r[doseKey], ler: r[lerKey] }));
  if (pts.length < 3) return { success: false, message: "Need at least 3 valid data points." };

  // Least squares fit: log(LER - B_guess) = log(A) - alpha * log(dose)
  // Simplified: try alpha = 0.5 (shot noise) and fit A, B
  let bestR2 = -Infinity, bestA = 1, bestB = 0, bestAlpha = 0.5;
  for (let alpha = 0.3; alpha <= 1.0; alpha += 0.05) {
    for (let B = 0; B <= pts[pts.length - 1].ler * 0.9; B += 0.05) {
      let sumXY = 0, sumX2 = 0, sumY = 0, sumX = 0, n = 0;
      for (const p of pts) {
        if (p.ler - B <= 0) continue;
        const x = Math.log(p.dose), y = Math.log(p.ler - B);
        sumXY += x * y; sumX2 += x * x; sumX += x; sumY += y; n++;
      }
      if (n < 2) continue;
      const slopeF = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
      const interceptF = (sumY - slopeF * sumX) / n;
      const A = Math.exp(interceptF);
      // R²
      const meanY = pts.reduce((s, p) => s + p.ler, 0) / pts.length;
      let ssTot = 0, ssRes = 0;
      for (const p of pts) {
        const pred = A * Math.pow(p.dose, -alpha) + B;
        ssRes += (p.ler - pred) ** 2;
        ssTot += (p.ler - meanY) ** 2;
      }
      const R2 = 1 - ssRes / Math.max(1e-10, ssTot);
      if (R2 > bestR2) { bestR2 = R2; bestA = A; bestB = B; bestAlpha = alpha; }
    }
  }

  // Generate fit curve
  const minD = Math.min(...pts.map(p => p.dose)), maxD = Math.max(...pts.map(p => p.dose));
  const fitCurve = [];
  for (let d = minD; d <= maxD * 1.3; d += (maxD - minD) / 50) {
    fitCurve.push({ dose: +d.toFixed(1), lerFit: +(bestA * Math.pow(d, -bestAlpha) + bestB).toFixed(3) });
  }

  // Merge experimental points
  const chartData = fitCurve.map(f => {
    const closest = pts.reduce((best, p) => Math.abs(p.dose - f.dose) < Math.abs(best.dose - f.dose) ? p : best, pts[0]);
    return { ...f, lerExp: Math.abs(closest.dose - f.dose) < (maxD - minD) / 30 ? closest.ler : undefined };
  });

  // Extract key metrics
  const doseForLER2 = bestA > 0 ? Math.pow((2 - bestB) / bestA, -1 / bestAlpha) : null;
  const doseForLER15 = bestA > 0 ? Math.pow((1.5 - bestB) / bestA, -1 / bestAlpha) : null;

  return {
    success: true, A: bestA, B: bestB, alpha: bestAlpha, R2: bestR2,
    chartData, pts, doseForLER2, doseForLER15,
    equation: `LER = ${bestA.toFixed(2)} × dose^(-${bestAlpha.toFixed(2)}) + ${bestB.toFixed(2)}`,
  };
}

// ════════════════════════════════════════════════════════════════
// UI PALETTE & STYLES
// ════════════════════════════════════════════════════════════════

const P = {
  bg: "#060b18", s1: "#0c1426", s2: "#131f38", bd: "#1e3058", bdL: "#2a4270",
  ac: "#00d4ff", acD: "#007799", euv: "#a855f7", sn: "#fbbf24", 
  warn: "#f59e0b", err: "#ef4444", ok: "#22c55e",
  t1: "#e2e8f0", t2: "#94a3b8", t3: "#475569",
};

const S = {
  root: { background: P.bg, minHeight: "100vh", color: P.t1, fontFamily: "'IBM Plex Mono', 'Cascadia Code', 'JetBrains Mono', monospace", fontSize: 12, lineHeight: 1.5 },
  hdr: { background: `linear-gradient(135deg, ${P.s1}, #080e22)`, borderBottom: `1px solid ${P.bd}`, padding: "16px 20px 12px", display: "flex", justifyContent: "space-between", alignItems: "flex-end" },
  title: { fontSize: 18, fontWeight: 800, margin: 0, background: `linear-gradient(90deg, ${P.ac}, ${P.euv})`, WebkitBackgroundClip: "text", WebkitTextFillColor: "transparent", letterSpacing: "0.3px" },
  sub: { fontSize: 9, color: P.t3, letterSpacing: "2px", textTransform: "uppercase", marginTop: 3 },
  grid: { display: "grid", gridTemplateColumns: "290px 1fr", minHeight: "calc(100vh - 70px)" },
  side: { background: P.s1, borderRight: `1px solid ${P.bd}`, padding: 14, overflowY: "auto", maxHeight: "calc(100vh - 70px)" },
  main: { padding: 14, overflowY: "auto", maxHeight: "calc(100vh - 70px)" },
  sec: { marginBottom: 12 },
  secT: { fontSize: 9, fontWeight: 700, textTransform: "uppercase", letterSpacing: "2px", color: P.ac, marginBottom: 8, display: "flex", alignItems: "center", gap: 5 },
  card: { background: P.s2, border: `1px solid ${P.bd}`, borderRadius: 5, padding: 12, marginBottom: 10 },
  inp: { width: "100%", background: P.s1, border: `1px solid ${P.bd}`, borderRadius: 3, padding: "5px 7px", color: P.t1, fontSize: 11, fontFamily: "inherit", outline: "none", boxSizing: "border-box" },
  sl: { width: "100%", accentColor: P.ac, margin: "1px 0" },
  mr: { display: "flex", justifyContent: "space-between", alignItems: "baseline", padding: "3px 0", borderBottom: `1px solid ${P.bd}15` },
  tab: (a) => ({ padding: "7px 12px", fontSize: 10, fontWeight: a ? 700 : 500, color: a ? P.ac : P.t2, background: a ? `${P.ac}11` : "transparent", border: "none", borderBottom: a ? `2px solid ${P.ac}` : "2px solid transparent", cursor: "pointer", fontFamily: "inherit", letterSpacing: "0.5px" }),
  alert: (t) => ({ padding: "8px 12px", borderRadius: 4, fontSize: 10, lineHeight: 1.5, marginBottom: 10, background: `${t === "w" ? P.warn : t === "e" ? P.err : P.ok}11`, border: `1px solid ${t === "w" ? P.warn : t === "e" ? P.err : P.ok}33`, color: t === "w" ? P.warn : t === "e" ? P.err : P.ok }),
  dropzone: { border: `2px dashed ${P.bd}`, borderRadius: 6, padding: "24px 16px", textAlign: "center", cursor: "pointer", transition: "all 0.2s", background: P.s1 },
  dropzoneActive: { border: `2px dashed ${P.ac}`, background: `${P.ac}08` },
  ttStyle: { background: P.s1, border: `1px solid ${P.bd}`, fontSize: 10, fontFamily: "inherit" },
};

function Sl({ label, unit, value, onChange, min, max, step }) {
  return (
    <div style={{ marginBottom: 6 }}>
      <div style={{ display: "flex", justifyContent: "space-between", fontSize: 10, color: P.t2, marginBottom: 2 }}>
        <span>{label}</span><span style={{ color: P.ac, fontWeight: 600 }}>{value}{unit && ` ${unit}`}</span>
      </div>
      <input type="range" min={min} max={max} step={step} value={value} onChange={e => onChange(+e.target.value)} style={S.sl} />
    </div>
  );
}

function Mv({ l, v, u, c, exp }) {
  const d = exp ? parseFloat(v).toExponential(2) : typeof v === "number" ? v.toFixed(3) : v;
  return (
    <div style={S.mr}>
      <span style={{ fontSize: 10, color: P.t2 }}>{l}</span>
      <span style={{ fontSize: 12, fontWeight: 600, color: c || P.t1 }}>{d} <span style={{ fontSize: 9, color: P.t3 }}>{u}</span></span>
    </div>
  );
}

function Badge({ color, children }) {
  return <span style={{ display: "inline-block", padding: "1px 7px", borderRadius: 3, fontSize: 9, fontWeight: 700, background: `${color}22`, color, border: `1px solid ${color}44` }}>{children}</span>;
}

// ════════════════════════════════════════════════════════════════
// MAIN APP
// ════════════════════════════════════════════════════════════════

export default function App() {
  const [nSn, setNSn] = useState(8);
  const [molSz, setMolSz] = useState(1.2);
  const [ringD, setRingD] = useState(1.5);
  const [filmTh, setFilmTh] = useState(12);
  const [dose, setDose] = useState(30);
  const [NA, setNA] = useState(0.55);
  const [hp, setHp] = useState(10);
  const [seB, setSeB] = useState(2);
  const [chemB, setChemB] = useState(1.5);
  const [thRE, setThRE] = useState(2);
  const [thSCF, setThSCF] = useState(14);
  const [agT, setAgT] = useState(24);
  const [temp, setTemp] = useState(23);
  const [hum, setHum] = useState(45);
  const [TON, setTON] = useState(1);
  const [tab, setTab] = useState("overview");

  // Data import states
  const [imgResult, setImgResult] = useState(null);
  const [csvResult, setCsvResult] = useState(null);
  const [imgScale, setImgScale] = useState(0.5); // nm/pixel
  const [dragOver, setDragOver] = useState(false);
  const canvasRef = useRef(null);
  const fileRef = useRef(null);
  const csvRef = useRef(null);

  const r = useMemo(() => {
    try { return runModel({ nSn, molSz, ringD, filmTh, dose, NA, hp, seB, chemB, thRE, thSCF, agT, temp, hum, TON }); }
    catch (e) { console.error(e); return null; }
  }, [nSn, molSz, ringD, filmTh, dose, NA, hp, seB, chemB, thRE, thSCF, agT, temp, hum, TON]);

  // Image processing handler
  const handleImage = useCallback((file) => {
    if (!file || !file.type.startsWith("image/")) return;
    const reader = new FileReader();
    reader.onload = (e) => {
      const img = new Image();
      img.onload = () => {
        const canvas = canvasRef.current;
        if (!canvas) return;
        const maxW = 512, scale = img.width > maxW ? maxW / img.width : 1;
        const w = Math.round(img.width * scale), h = Math.round(img.height * scale);
        canvas.width = w; canvas.height = h;
        const ctx = canvas.getContext("2d");
        ctx.drawImage(img, 0, 0, w, h);
        const imgData = ctx.getImageData(0, 0, w, h).data;
        const result = analyzeImageForLER(imgData, w, h);
        if (result.success) {
          result.lerNm = result.lerAvg * imgScale;
          result.lwrNm = result.lwr * imgScale;
          result.cdNm = result.cdMean * imgScale;
          result.fileName = file.name;
        }
        setImgResult(result);
      };
      img.src = e.target.result;
    };
    reader.readAsDataURL(file);
  }, [imgScale]);

  const handleCSV = useCallback((file) => {
    const reader = new FileReader();
    reader.onload = (e) => {
      const parsed = parseCSV(e.target.result);
      if (parsed) {
        const fit = fitDoseLER(parsed);
        setCsvResult(fit);
      }
    };
    reader.readAsText(file);
  }, []);

  const onDrop = useCallback((e) => {
    e.preventDefault(); setDragOver(false);
    const file = e.dataTransfer?.files?.[0];
    if (!file) return;
    if (file.type.startsWith("image/")) handleImage(file);
    else if (file.name.endsWith(".csv") || file.name.endsWith(".tsv") || file.name.endsWith(".txt")) handleCSV(file);
  }, [handleImage, handleCSV]);

  if (!r) return <div style={S.root}><p style={{ padding: 20 }}>Calculation error.</p></div>;

  const tabs = [
    ["overview", "Overview"], ["aging", "Aging"], ["dose", "Dose"], ["profile", "DR Profile"],
    ["window", "Process Window"], ["data", "Data Import"], ["deploy", "Deploy Guide"]
  ];

  return (
    <div style={S.root}>
      <div style={S.hdr}>
        <div>
          <h1 style={S.title}>EUV MOR-PR Stochastic Aging Simulator v2</h1>
          <div style={S.sub}>Sn nano-ring MOR · Fukuda network model · experimental data integration · DFT-free</div>
        </div>
        <div style={{ display: "flex", gap: 6, alignItems: "center" }}>
          <Badge color={P.ok}>LIVE</Badge>
          <Badge color={P.euv}>13.5 nm</Badge>
        </div>
      </div>

      <div style={S.grid}>
        {/* SIDEBAR */}
        <div style={S.side}>
          <div style={S.sec}>
            <div style={S.secT}><span style={{ color: P.sn }}>◆</span> Molecular Structure</div>
            <Sl label="Sn atoms/unit" value={nSn} onChange={setNSn} min={2} max={20} step={1} />
            <Sl label="Mol. size" unit="nm" value={molSz} onChange={setMolSz} min={0.5} max={3} step={0.1} />
            <Sl label="Ring ∅" unit="nm" value={ringD} onChange={setRingD} min={0.8} max={4} step={0.1} />
          </div>
          <div style={S.sec}>
            <div style={S.secT}><span style={{ color: P.euv }}>◆</span> Process</div>
            <Sl label="Film" unit="nm" value={filmTh} onChange={setFilmTh} min={5} max={40} step={1} />
            <Sl label="Dose" unit="mJ/cm²" value={dose} onChange={setDose} min={5} max={120} step={1} />
            <Sl label="NA" value={NA} onChange={setNA} min={0.25} max={0.75} step={0.01} />
            <Sl label="Half-pitch" unit="nm" value={hp} onChange={setHp} min={5} max={30} step={0.5} />
            <Sl label="TON" value={TON} onChange={setTON} min={1} max={10} step={1} />
          </div>
          <div style={S.sec}>
            <div style={S.secT}><span style={{ color: P.ac }}>◆</span> Stochastics</div>
            <Sl label="SE blur" unit="nm" value={seB} onChange={setSeB} min={0.5} max={8} step={0.5} />
            <Sl label="Chem blur" unit="nm" value={chemB} onChange={setChemB} min={0.5} max={8} step={0.5} />
            <Sl label="th_RE" value={thRE} onChange={setThRE} min={1} max={5} step={1} />
            <Sl label="th_SCF" unit="/27" value={thSCF} onChange={setThSCF} min={5} max={22} step={1} />
          </div>
          <div style={S.sec}>
            <div style={S.secT}><span style={{ color: P.warn }}>◆</span> Aging</div>
            <Sl label="Queue" unit="h" value={agT} onChange={setAgT} min={0} max={720} step={1} />
            <Sl label="Temp" unit="°C" value={temp} onChange={setTemp} min={10} max={50} step={1} />
            <Sl label="RH" unit="%" value={hum} onChange={setHum} min={10} max={90} step={5} />
          </div>
          {/* Quick summary box */}
          <div style={{ ...S.card, background: `${P.ac}08`, border: `1px solid ${P.ac}22` }}>
            <div style={{ fontSize: 9, color: P.ac, fontWeight: 700, marginBottom: 4, letterSpacing: "1px" }}>QUICK SUMMARY</div>
            <div style={{ fontSize: 11, color: P.t2, lineHeight: 1.6 }}>
              k₁ = <b style={{ color: r.k1 < 0.3 ? P.err : P.t1 }}>{r.k1.toFixed(3)}</b> · LER = <b style={{ color: r.ler > 2.5 ? P.err : P.t1 }}>{r.ler.toFixed(2)}</b> nm<br />
              Aged LER = <b style={{ color: r.agLER > r.ler * 1.2 ? P.warn : P.t1 }}>{r.agLER.toFixed(2)}</b> nm · Dose ×{r.doseCor.toFixed(2)}
            </div>
          </div>
        </div>

        {/* MAIN */}
        <div style={S.main}>
          {r.ler > 2.5 && <div style={S.alert("e")}>⚠ LER ({r.ler.toFixed(2)} nm) exceeds 2.5 nm spec. Increase dose or thresholds.</div>}
          {r.k1 < 0.3 && <div style={S.alert("w")}>⚡ k₁ = {r.k1.toFixed(3)} — aggressive regime. Stochastic defects may spike.</div>}
          {r.ag.intact < 0.9 && <div style={S.alert("w")}>⏳ {((1 - r.ag.intact) * 100).toFixed(1)}% aging @ {agT}h. Dose correction ×{r.doseCor.toFixed(2)}</div>}

          <div style={{ display: "flex", gap: 2, marginBottom: 14, borderBottom: `1px solid ${P.bd}`, flexWrap: "wrap" }}>
            {tabs.map(([k, l]) => <button key={k} style={S.tab(tab === k)} onClick={() => setTab(k)}>{l}</button>)}
          </div>

          <canvas ref={canvasRef} style={{ display: "none" }} />

          {/* ═══ OVERVIEW ═══ */}
          {tab === "overview" && (
            <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 10 }}>
              <div style={S.card}>
                <div style={{ ...S.secT, marginBottom: 6 }}>Lithography (Fresh)</div>
                <Mv l="k₁" v={r.k1} c={r.k1 < 0.3 ? P.err : r.k1 < 0.35 ? P.warn : P.ok} />
                <Mv l="LER (3σ)" v={r.ler} u="nm" c={r.ler > 2.5 ? P.err : r.ler > 1.5 ? P.warn : P.ok} />
                <Mv l="LWR (3σ)" v={r.lwr} u="nm" />
                <Mv l="Univ. DR log-slope" v={r.uLS} />
                <Mv l="Image contrast" v={r.imgC} />
                <Mv l="Defect prob" v={r.defP} exp c={r.defP > 0.01 ? P.err : P.ok} />
              </div>
              <div style={S.card}>
                <div style={{ ...S.secT, marginBottom: 6 }}>Photochemistry</div>
                <Mv l="μ (abs. coeff)" v={r.mu} u="nm⁻¹" c={P.euv} />
                <Mv l="Photons/mol" v={r.ppm} />
                <Mv l="Reactions/mol" v={r.rpm} />
                <Mv l="p(MSS) exposed" v={r.pMSS_e} c={P.ok} />
                <Mv l="p(MSS) unexposed" v={r.pMSS_u} c={r.pMSS_u > 0.1 ? P.err : P.ok} />
                <Mv l="p(SC) exposed" v={r.pSC_e} c={P.ok} />
              </div>
              <div style={S.card}>
                <div style={{ ...S.secT, marginBottom: 6 }}>Aging @ {agT}h</div>
                <Mv l="Intact %" v={r.ag.intact * 100} u="%" c={r.ag.intact < 0.9 ? P.warn : P.ok} />
                <Mv l="Aged LER" v={r.agLER} u="nm" c={r.agLER > r.ler * 1.2 ? P.warn : P.t1} />
                <Mv l="Aged LWR" v={r.agLWR} u="nm" />
                <Mv l="Dose correction" v={r.doseCor} u="×" c={r.doseCor > 1.1 ? P.warn : P.t1} />
                <Mv l="Overdispersion Δ" v={r.ag.odIncrease} />
              </div>
              <div style={S.card}>
                <div style={{ ...S.secT, marginBottom: 6 }}>Process Score</div>
                <ResponsiveContainer width="100%" height={200}>
                  <RadarChart data={r.radar}>
                    <PolarGrid stroke={P.bd} />
                    <PolarAngleAxis dataKey="m" tick={{ fill: P.t2, fontSize: 8 }} />
                    <PolarRadiusAxis angle={90} domain={[0, 100]} tick={false} axisLine={false} />
                    <Radar dataKey="v" stroke={P.ac} fill={P.ac} fillOpacity={0.2} strokeWidth={2} />
                  </RadarChart>
                </ResponsiveContainer>
              </div>
            </div>
          )}

          {/* ═══ AGING ═══ */}
          {tab === "aging" && (
            <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 10 }}>
              <div style={S.card}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>LER / LWR vs Queue Time</div>
                <ResponsiveContainer width="100%" height={230}>
                  <LineChart data={r.agTL}>
                    <CartesianGrid strokeDasharray="3 3" stroke={P.bd} />
                    <XAxis dataKey="tLabel" tick={{ fill: P.t2, fontSize: 9 }} />
                    <YAxis tick={{ fill: P.t2, fontSize: 9 }} />
                    <Tooltip contentStyle={S.ttStyle} />
                    <Line type="monotone" dataKey="ler" stroke={P.ac} name="LER" strokeWidth={2} dot={false} />
                    <Line type="monotone" dataKey="lwr" stroke={P.euv} name="LWR" strokeWidth={2} dot={false} />
                    <ReferenceLine y={2.5} stroke={P.err} strokeDasharray="4 4" />
                    <Legend wrapperStyle={{ fontSize: 9 }} />
                  </LineChart>
                </ResponsiveContainer>
              </div>
              <div style={S.card}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>Molecular Integrity</div>
                <ResponsiveContainer width="100%" height={230}>
                  <ComposedChart data={r.agTL}>
                    <CartesianGrid strokeDasharray="3 3" stroke={P.bd} />
                    <XAxis dataKey="tLabel" tick={{ fill: P.t2, fontSize: 9 }} />
                    <YAxis tick={{ fill: P.t2, fontSize: 9 }} />
                    <Tooltip contentStyle={S.ttStyle} />
                    <Area type="monotone" dataKey="intact" stroke={P.ok} fill={P.ok} fillOpacity={0.1} name="Intact %" strokeWidth={2} />
                    <Line type="monotone" dataKey="absPct" stroke={P.sn} name="Abs %" strokeWidth={2} dot={false} />
                    <Legend wrapperStyle={{ fontSize: 9 }} />
                  </ComposedChart>
                </ResponsiveContainer>
              </div>
              <div style={S.card}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>Dose Correction Factor</div>
                <ResponsiveContainer width="100%" height={230}>
                  <AreaChart data={r.agTL}>
                    <CartesianGrid strokeDasharray="3 3" stroke={P.bd} />
                    <XAxis dataKey="tLabel" tick={{ fill: P.t2, fontSize: 9 }} />
                    <YAxis tick={{ fill: P.t2, fontSize: 9 }} />
                    <Tooltip contentStyle={S.ttStyle} />
                    <Area type="monotone" dataKey="doseCor" stroke={P.warn} fill={P.warn} fillOpacity={0.1} name="×factor" strokeWidth={2} />
                  </AreaChart>
                </ResponsiveContainer>
              </div>
              <div style={S.card}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>Aging Recommendation</div>
                <div style={{ fontSize: 10, lineHeight: 1.7, color: P.t2 }}>
                  {r.ag.intact > 0.95 ? <span style={{ color: P.ok }}>✓ Queue {agT}h @ {temp}°C/{hum}%RH — minimal degradation.</span>
                    : r.ag.intact > 0.85 ? <span style={{ color: P.warn }}>⚠ Moderate aging ({((1 - r.ag.intact) * 100).toFixed(1)}%). Apply dose ×{r.doseCor.toFixed(2)}.</span>
                      : <span style={{ color: P.err }}>✗ Severe aging. Reduce queue to {"<"}{Math.round(agT * 0.4)}h or cool to {"<"}18°C.</span>}
                  <br />Max queue (95% intact): ~{(() => {
                    for (const t of [1, 2, 4, 8, 12, 24, 48, 72, 96, 120, 168, 336, 720]) {
                      if (agingModel(t, temp, hum, {}).intact < 0.95) return t < 24 ? `${t}h` : `${Math.round(t / 24)}d`;
                    } return ">30d";
                  })()}
                </div>
              </div>
            </div>
          )}

          {/* ═══ DOSE ═══ */}
          {tab === "dose" && (
            <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 10 }}>
              <div style={S.card}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>LER vs Dose</div>
                <ResponsiveContainer width="100%" height={260}>
                  <LineChart data={r.dSweep}>
                    <CartesianGrid strokeDasharray="3 3" stroke={P.bd} />
                    <XAxis dataKey="dose" tick={{ fill: P.t2, fontSize: 9 }} />
                    <YAxis tick={{ fill: P.t2, fontSize: 9 }} />
                    <Tooltip contentStyle={S.ttStyle} />
                    <Line type="monotone" dataKey="ler" stroke={P.ac} name="LER" strokeWidth={2} dot={false} />
                    <Line type="monotone" dataKey="lwr" stroke={P.euv} name="LWR" strokeWidth={2} dot={false} />
                    <ReferenceLine y={2.5} stroke={P.err} strokeDasharray="4 4" />
                    <ReferenceLine x={dose} stroke={P.sn} strokeDasharray="3 3" />
                    <Legend wrapperStyle={{ fontSize: 9 }} />
                  </LineChart>
                </ResponsiveContainer>
              </div>
              <div style={S.card}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>Contrast vs Dose</div>
                <ResponsiveContainer width="100%" height={260}>
                  <AreaChart data={r.dSweep}>
                    <CartesianGrid strokeDasharray="3 3" stroke={P.bd} />
                    <XAxis dataKey="dose" tick={{ fill: P.t2, fontSize: 9 }} />
                    <YAxis tick={{ fill: P.t2, fontSize: 9 }} domain={[0, 1]} />
                    <Tooltip contentStyle={S.ttStyle} />
                    <Area type="monotone" dataKey="contrast" stroke={P.ok} fill={P.ok} fillOpacity={0.15} strokeWidth={2} />
                    <ReferenceLine x={dose} stroke={P.sn} strokeDasharray="3 3" />
                  </AreaChart>
                </ResponsiveContainer>
              </div>
            </div>
          )}

          {/* ═══ PROFILE ═══ */}
          {tab === "profile" && (
            <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 10 }}>
              <div style={S.card}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>Dissolution Rate Profile</div>
                <ResponsiveContainer width="100%" height={260}>
                  <ComposedChart data={r.drProf}>
                    <CartesianGrid strokeDasharray="3 3" stroke={P.bd} />
                    <XAxis dataKey="x" tick={{ fill: P.t2, fontSize: 9 }} />
                    <YAxis tick={{ fill: P.t2, fontSize: 9 }} />
                    <Tooltip contentStyle={S.ttStyle} />
                    <Area type="monotone" dataKey="dr" stroke={P.ac} fill={P.ac} fillOpacity={0.12} name="DR" strokeWidth={2} />
                    <Line type="monotone" dataKey="drV" stroke={P.warn} name="Variance" strokeWidth={1.5} dot={false} strokeDasharray="4 2" />
                    <Legend wrapperStyle={{ fontSize: 9 }} />
                  </ComposedChart>
                </ResponsiveContainer>
              </div>
              <div style={S.card}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>Image Intensity</div>
                <ResponsiveContainer width="100%" height={260}>
                  <AreaChart data={r.drProf}>
                    <CartesianGrid strokeDasharray="3 3" stroke={P.bd} />
                    <XAxis dataKey="x" tick={{ fill: P.t2, fontSize: 9 }} />
                    <YAxis tick={{ fill: P.t2, fontSize: 9 }} domain={[0, 1]} />
                    <Tooltip contentStyle={S.ttStyle} />
                    <Area type="monotone" dataKey="img" stroke={P.euv} fill={P.euv} fillOpacity={0.12} strokeWidth={2} />
                    <ReferenceLine y={0.5} stroke={P.t3} strokeDasharray="4 4" />
                  </AreaChart>
                </ResponsiveContainer>
              </div>
              <div style={{ ...S.card, gridColumn: "1/-1" }}>
                <div style={{ display: "grid", gridTemplateColumns: "repeat(4,1fr)", gap: 12 }}>
                  {[
                    ["DR SLOPE", r.drSlope.toFixed(3), "nm⁻¹"],
                    ["UNIV. LOG-SLOPE", r.uLS.toFixed(1), ""],
                    ["MOL LAYERS", r.nMolTh, ""],
                    ["TOTAL BLUR", r.blur.toFixed(1), "nm"],
                  ].map(([l, v, u], i) => (
                    <div key={i}>
                      <div style={{ fontSize: 9, color: P.t3, marginBottom: 3 }}>{l}</div>
                      <div style={{ fontSize: 16, fontWeight: 700, color: P.ac }}>{v}</div>
                      <div style={{ fontSize: 9, color: P.t3 }}>{u}</div>
                    </div>
                  ))}
                </div>
              </div>
            </div>
          )}

          {/* ═══ PROCESS WINDOW ═══ */}
          {tab === "window" && (
            <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 10 }}>
              <div style={S.card}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>Radar</div>
                <ResponsiveContainer width="100%" height={260}>
                  <RadarChart data={r.radar}>
                    <PolarGrid stroke={P.bd} />
                    <PolarAngleAxis dataKey="m" tick={{ fill: P.t2, fontSize: 9 }} />
                    <PolarRadiusAxis angle={90} domain={[0, 100]} tick={{ fill: P.t3, fontSize: 8 }} />
                    <Radar dataKey="v" stroke={P.ac} fill={P.ac} fillOpacity={0.2} strokeWidth={2} />
                  </RadarChart>
                </ResponsiveContainer>
              </div>
              <div style={S.card}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 8 }}>Score Breakdown</div>
                {r.radar.map((d, i) => (
                  <div key={i} style={{ display: "flex", alignItems: "center", gap: 6, marginBottom: 5 }}>
                    <div style={{ width: 90, fontSize: 10, color: P.t2 }}>{d.m}</div>
                    <div style={{ flex: 1, height: 6, background: P.bd, borderRadius: 3, overflow: "hidden" }}>
                      <div style={{ width: `${d.v}%`, height: "100%", borderRadius: 3, background: d.v > 70 ? P.ok : d.v > 40 ? P.warn : P.err }} />
                    </div>
                    <div style={{ width: 30, fontSize: 10, fontWeight: 700, textAlign: "right", color: d.v > 70 ? P.ok : d.v > 40 ? P.warn : P.err }}>{d.v}%</div>
                  </div>
                ))}
                <div style={{ marginTop: 12, padding: 10, background: `${P.ac}08`, borderRadius: 4, border: `1px solid ${P.ac}22`, textAlign: "center" }}>
                  <div style={{ fontSize: 9, color: P.ac, fontWeight: 700 }}>COMPOSITE</div>
                  <div style={{ fontSize: 22, fontWeight: 800, color: P.ac }}>{(r.radar.reduce((s, d) => s + (+d.v), 0) / r.radar.length).toFixed(0)}<span style={{ fontSize: 11, color: P.t3 }}>/100</span></div>
                </div>
              </div>
              <div style={{ ...S.card, gridColumn: "1/-1" }}>
                <div style={S.secT}>Action Items</div>
                <div style={{ fontSize: 10, lineHeight: 1.8, color: P.t2 }}>
                  {r.k1 < 0.35 && <div>▸ <span style={{ color: P.warn }}>Resolution:</span> k₁={r.k1.toFixed(3)}. Apply OPC/SMO.</div>}
                  {r.ler > 2.0 && <div>▸ <span style={{ color: P.warn }}>LER:</span> {r.ler.toFixed(2)} nm → increase dose to ~{Math.round(dose * (r.ler / 1.5) ** 2)} mJ/cm².</div>}
                  {r.defP > 1e-4 && <div>▸ <span style={{ color: P.err }}>Defects:</span> P={r.defP.toExponential(1)}. Tighten thresholds.</div>}
                  {r.ag.intact < 0.95 && <div>▸ <span style={{ color: P.warn }}>Aging:</span> enforce FIFO, target {"<"}{Math.round(agT * 0.5)}h queue.</div>}
                  <div>▸ <span style={{ color: P.ac }}>Blur:</span> {r.blur.toFixed(1)} nm / {hp} nm HP = {(r.blur / hp * 100).toFixed(0)}%. Target {"<"}30%.</div>
                </div>
              </div>
            </div>
          )}

          {/* ═══ DATA IMPORT ═══ */}
          {tab === "data" && (
            <div>
              <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 10 }}>
                {/* Image Upload */}
                <div style={S.card}>
                  <div style={S.secT}>SEM / TEM / AFM Image Analysis</div>
                  <div style={{ fontSize: 10, color: P.t2, marginBottom: 8, lineHeight: 1.6 }}>
                    Drop or select a line/space pattern SEM image. The system will auto-detect edges via Otsu thresholding, measure LER/LWR in pixels, and convert using your scale factor.
                  </div>
                  <div style={{ marginBottom: 8 }}>
                    <Sl label="Image scale" unit="nm/px" value={imgScale} onChange={setImgScale} min={0.1} max={5} step={0.1} />
                  </div>
                  <div
                    style={{ ...S.dropzone, ...(dragOver ? S.dropzoneActive : {}) }}
                    onDragOver={e => { e.preventDefault(); setDragOver(true); }}
                    onDragLeave={() => setDragOver(false)}
                    onDrop={onDrop}
                    onClick={() => fileRef.current?.click()}
                  >
                    <div style={{ fontSize: 20, marginBottom: 6 }}>📷</div>
                    <div style={{ fontSize: 11, color: P.t2 }}>Drop SEM/TEM/AFM image here</div>
                    <div style={{ fontSize: 9, color: P.t3, marginTop: 4 }}>PNG, JPG, BMP · Line/Space patterns</div>
                  </div>
                  <input ref={fileRef} type="file" accept="image/*" style={{ display: "none" }}
                    onChange={e => e.target.files?.[0] && handleImage(e.target.files[0])} />

                  {imgResult && (
                    <div style={{ marginTop: 10 }}>
                      {imgResult.success ? (
                        <>
                          <div style={S.alert("ok")}>✓ Edge detection complete — {imgResult.nEdges} edges found</div>
                          <Mv l="LER (3σ, pixel)" v={imgResult.lerAvg.toFixed(1)} u="px" />
                          <Mv l="LER (3σ, scaled)" v={imgResult.lerNm.toFixed(2)} u="nm" c={imgResult.lerNm > 2.5 ? P.err : P.ok} />
                          <Mv l="LWR (3σ, scaled)" v={imgResult.lwrNm.toFixed(2)} u="nm" />
                          <Mv l="CD mean" v={imgResult.cdNm.toFixed(1)} u="nm" />
                          <Mv l="CD σ" v={(imgResult.cdStd * imgScale).toFixed(2)} u="nm" />
                          <Mv l="Otsu threshold" v={imgResult.threshold.toFixed(3)} />
                          <div style={{ fontSize: 10, color: P.t2, marginTop: 8 }}>
                            <b>Model vs Exp:</b> Model LER = {r.ler.toFixed(2)} nm, Exp LER = {imgResult.lerNm.toFixed(2)} nm
                            {Math.abs(r.ler - imgResult.lerNm) / imgResult.lerNm < 0.2
                              ? <Badge color={P.ok}> {"<"}20% deviation</Badge>
                              : <Badge color={P.warn}> {((Math.abs(r.ler - imgResult.lerNm) / imgResult.lerNm) * 100).toFixed(0)}% deviation</Badge>}
                          </div>
                        </>
                      ) : (
                        <div style={S.alert("e")}>{imgResult.message}</div>
                      )}
                    </div>
                  )}

                  {imgResult?.success && imgResult.edgeProfile && (
                    <div style={{ marginTop: 10 }}>
                      <div style={{ fontSize: 10, fontWeight: 700, marginBottom: 4 }}>Edge Profile (left edge)</div>
                      <ResponsiveContainer width="100%" height={150}>
                        <LineChart data={imgResult.edgeProfile.slice(0, 150)}>
                          <XAxis dataKey="y" tick={{ fill: P.t3, fontSize: 8 }} />
                          <YAxis tick={{ fill: P.t3, fontSize: 8 }} domain={['dataMin - 5', 'dataMax + 5']} />
                          <Line type="linear" dataKey="x" stroke={P.ac} strokeWidth={1} dot={false} name="Left edge" />
                          <Line type="linear" dataKey="xR" stroke={P.euv} strokeWidth={1} dot={false} name="Right edge" />
                          <Tooltip contentStyle={S.ttStyle} />
                        </LineChart>
                      </ResponsiveContainer>
                    </div>
                  )}

                  {imgResult?.success && imgResult.psd && (
                    <div style={{ marginTop: 8 }}>
                      <div style={{ fontSize: 10, fontWeight: 700, marginBottom: 4 }}>LER Power Spectrum</div>
                      <ResponsiveContainer width="100%" height={120}>
                        <LineChart data={imgResult.psd}>
                          <XAxis dataKey="freq" tick={{ fill: P.t3, fontSize: 8 }} />
                          <YAxis tick={{ fill: P.t3, fontSize: 8 }} scale="log" domain={['auto', 'auto']} />
                          <Line type="monotone" dataKey="power" stroke={P.sn} strokeWidth={1.5} dot={false} />
                          <Tooltip contentStyle={S.ttStyle} />
                        </LineChart>
                      </ResponsiveContainer>
                    </div>
                  )}
                </div>

                {/* CSV Upload */}
                <div style={S.card}>
                  <div style={S.secT}>Dose-LER / Spectroscopy CSV</div>
                  <div style={{ fontSize: 10, color: P.t2, marginBottom: 8, lineHeight: 1.6 }}>
                    Upload a CSV with columns like <code style={{ color: P.ac }}>dose, ler</code> (or <code style={{ color: P.ac }}>energy, roughness</code>).
                    The system fits LER = A·dose<sup>−α</sup> + B and extracts optimal dose targets.
                  </div>
                  <div style={{ fontSize: 9, color: P.t3, marginBottom: 8, padding: 8, background: P.s1, borderRadius: 3, fontFamily: "monospace" }}>
                    Example CSV:<br />
                    dose,ler<br />
                    10,3.8<br />
                    15,3.1<br />
                    20,2.6<br />
                    30,2.1<br />
                    50,1.6<br />
                    80,1.3
                  </div>
                  <div
                    style={{ ...S.dropzone, ...(dragOver ? S.dropzoneActive : {}) }}
                    onClick={() => csvRef.current?.click()}
                    onDragOver={e => { e.preventDefault(); setDragOver(true); }}
                    onDragLeave={() => setDragOver(false)}
                    onDrop={onDrop}
                  >
                    <div style={{ fontSize: 20, marginBottom: 6 }}>📊</div>
                    <div style={{ fontSize: 11, color: P.t2 }}>Drop CSV/TSV file here</div>
                  </div>
                  <input ref={csvRef} type="file" accept=".csv,.tsv,.txt" style={{ display: "none" }}
                    onChange={e => e.target.files?.[0] && handleCSV(e.target.files[0])} />

                  {csvResult && (
                    <div style={{ marginTop: 10 }}>
                      {csvResult.success ? (
                        <>
                          <div style={S.alert("ok")}>✓ Fit complete — R² = {csvResult.R2.toFixed(4)}</div>
                          <div style={{ fontSize: 10, color: P.t2, marginBottom: 6, fontFamily: "monospace" }}>{csvResult.equation}</div>
                          <Mv l="α (dose exponent)" v={csvResult.alpha.toFixed(3)} c={P.euv} />
                          <Mv l="A (coefficient)" v={csvResult.A.toFixed(3)} />
                          <Mv l="B (floor)" v={csvResult.B.toFixed(3)} u="nm" />
                          <Mv l="R²" v={csvResult.R2.toFixed(4)} c={csvResult.R2 > 0.95 ? P.ok : P.warn} />
                          {csvResult.doseForLER2 && <Mv l="Dose for LER<2.0" v={csvResult.doseForLER2.toFixed(1)} u="mJ/cm²" c={P.ac} />}
                          {csvResult.doseForLER15 && <Mv l="Dose for LER<1.5" v={csvResult.doseForLER15.toFixed(1)} u="mJ/cm²" c={P.ok} />}

                          <div style={{ marginTop: 10 }}>
                            <div style={{ fontSize: 10, fontWeight: 700, marginBottom: 4 }}>Fit vs Experimental</div>
                            <ResponsiveContainer width="100%" height={200}>
                              <ComposedChart data={csvResult.chartData}>
                                <CartesianGrid strokeDasharray="3 3" stroke={P.bd} />
                                <XAxis dataKey="dose" tick={{ fill: P.t2, fontSize: 9 }} />
                                <YAxis tick={{ fill: P.t2, fontSize: 9 }} />
                                <Tooltip contentStyle={S.ttStyle} />
                                <Line type="monotone" dataKey="lerFit" stroke={P.ac} name="Fit" strokeWidth={2} dot={false} />
                                <Scatter dataKey="lerExp" fill={P.sn} name="Exp" />
                                <ReferenceLine y={2.0} stroke={P.err} strokeDasharray="4 4" />
                                <Legend wrapperStyle={{ fontSize: 9 }} />
                              </ComposedChart>
                            </ResponsiveContainer>
                          </div>

                          <div style={{ fontSize: 10, color: P.t2, marginTop: 8, lineHeight: 1.6 }}>
                            {csvResult.alpha > 0.45 && csvResult.alpha < 0.55
                              ? <span style={{ color: P.ok }}>α ≈ 0.5 confirms photon shot noise dominance (Fukuda/Gallatin model).</span>
                              : csvResult.alpha > 0.55
                                ? <span style={{ color: P.warn }}>α {">"} 0.5 suggests additional stochastic sources beyond shot noise (SE clustering, material inhomogeneity).</span>
                                : <span style={{ color: P.warn }}>α {"<"} 0.45 may indicate dose-dependent chemical amplification or saturation effects.</span>}
                          </div>
                        </>
                      ) : (
                        <div style={S.alert("e")}>{csvResult.message}</div>
                      )}
                    </div>
                  )}
                </div>
              </div>
            </div>
          )}

          {/* ═══ DEPLOY GUIDE ═══ */}
          {tab === "deploy" && (
            <div style={{ maxWidth: 700 }}>
              <div style={S.card}>
                <div style={S.secT}>🚀 Deployment Guide</div>
                <div style={{ fontSize: 11, lineHeight: 2, color: P.t2 }}>

                  <div style={{ padding: "10px 12px", background: P.s1, borderRadius: 4, marginBottom: 12, border: `1px solid ${P.bd}` }}>
                    <div style={{ color: P.ac, fontWeight: 700, marginBottom: 6, fontSize: 12 }}>Option 1: Vercel (추천 — 가장 쉬움)</div>
                    <div>1. <a href="https://github.com" style={{ color: P.ac }}>github.com</a>에 새 repo 생성</div>
                    <div>2. <code style={{ color: P.sn, background: `${P.sn}11`, padding: "1px 4px", borderRadius: 2 }}>npx create-react-app euv-simulator</code> 로 프로젝트 생성</div>
                    <div>3. 이 코드를 <code style={{ color: P.sn }}>src/App.jsx</code>에 붙여넣기</div>
                    <div>4. <code style={{ color: P.sn }}>npm install recharts mathjs lodash</code></div>
                    <div>5. GitHub에 push</div>
                    <div>6. <a href="https://vercel.com" style={{ color: P.ac }}>vercel.com</a> → Import → GitHub repo 선택 → Deploy</div>
                    <div>7. 자동으로 <code style={{ color: P.ok }}>https://euv-simulator.vercel.app</code> 같은 URL 발급</div>
                  </div>

                  <div style={{ padding: "10px 12px", background: P.s1, borderRadius: 4, marginBottom: 12, border: `1px solid ${P.bd}` }}>
                    <div style={{ color: P.euv, fontWeight: 700, marginBottom: 6, fontSize: 12 }}>Option 2: Netlify</div>
                    <div>Vercel과 동일한 절차. <a href="https://netlify.com" style={{ color: P.ac }}>netlify.com</a>에서 GitHub repo 연결.</div>
                    <div>무료 플랜으로 월 100GB 트래픽, 커스텀 도메인 지원.</div>
                  </div>

                  <div style={{ padding: "10px 12px", background: P.s1, borderRadius: 4, marginBottom: 12, border: `1px solid ${P.bd}` }}>
                    <div style={{ color: P.sn, fontWeight: 700, marginBottom: 6, fontSize: 12 }}>Option 3: 사내 서버</div>
                    <div>1. <code style={{ color: P.sn }}>npm run build</code> 로 정적 파일 생성</div>
                    <div>2. <code style={{ color: P.sn }}>build/</code> 폴더를 사내 웹서버 (Apache/Nginx)에 복사</div>
                    <div>3. 사내 인트라넷 URL로 접근 가능</div>
                    <div>4. 보안이 필요하면 사내 VPN 뒤에 배치</div>
                  </div>

                  <div style={{ padding: "10px 12px", background: P.s1, borderRadius: 4, border: `1px solid ${P.bd}` }}>
                    <div style={{ color: P.ok, fontWeight: 700, marginBottom: 6, fontSize: 12 }}>Option 4: 단일 HTML (가장 단순)</div>
                    <div>이 앱을 단일 HTML 파일로 변환하면 어떤 웹서버에든 올릴 수 있습니다.</div>
                    <div>USB에 담아서 전달하거나, 이메일 첨부도 가능합니다.</div>
                    <div style={{ marginTop: 4, color: P.t3 }}>필요시 요청해 주세요 — HTML 버전으로 변환해 드립니다.</div>
                  </div>

                </div>
              </div>

              <div style={S.card}>
                <div style={S.secT}>📦 package.json dependencies</div>
                <pre style={{ fontSize: 10, color: P.sn, background: P.s1, padding: 10, borderRadius: 4, overflow: "auto" }}>
{`{
  "dependencies": {
    "react": "^18.2.0",
    "react-dom": "^18.2.0",
    "recharts": "^2.12.0",
    "mathjs": "^12.0.0",
    "lodash": "^4.17.21"
  }
}`}
                </pre>
              </div>

              <div style={S.card}>
                <div style={S.secT}>🔬 Physics Model References</div>
                <div style={{ fontSize: 10, color: P.t2, lineHeight: 1.8 }}>
                  <div>• Fukuda, J. Appl. Phys. 137, 204902 (2025) — Directional network model</div>
                  <div>• Binomial → Beta-binomial → Mixed Bernoulli convolution PMF hierarchy</div>
                  <div>• Arrhenius aging kinetics for Sn-O-Sn hydrolysis + aggregation</div>
                  <div>• LER = ΔDR / (∂DR/∂x) with k₁-dependent correction</div>
                  <div>• Stochastic defect prob = ∫ pmf(DR,r) dDR from threshold to ∞</div>
                  <div>• Universal log-slope = (1/Y)(∂Y/∂x)(λ/NA)</div>
                </div>
              </div>
            </div>
          )}

        </div>
      </div>
    </div>
  );
}
