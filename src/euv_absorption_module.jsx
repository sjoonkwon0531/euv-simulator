import { useState, useMemo, useCallback, useRef } from "react";
import * as math from "mathjs";
import {
  LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend,
  ResponsiveContainer, AreaChart, Area, ComposedChart, Bar, BarChart,
  ScatterChart, Scatter, ReferenceLine
} from "recharts";

// ════════════════════════════════════════════════════════════════
// CXRO / Henke Atomic Scattering Factor Database (f₂ values)
// Source: Henke, Gullikson & Davis, At. Data Nucl. Data Tables 54(2), 181-342 (1993)
// Updated values from https://henke.lbl.gov/optical_constants/
//
// f₂ relates to the photoabsorption cross-section:
//   σ_a(E) = 2 · r_e · λ · f₂(E)
// where r_e = 2.818 × 10⁻¹³ cm (classical electron radius)
//
// We tabulate f₂ for each element at energies spanning 80-105 eV (λ ≈ 12-15.5 nm)
// ════════════════════════════════════════════════════════════════

const r_e = 2.8179e-13; // classical electron radius (cm)
const hc = 1.23984e-4;  // h·c in eV·cm → hc/E(eV) = λ(cm)
const N_A = 6.022e23;

// f₂ values from Henke tables for key elements at EUV energies
// Format: { E(eV): f₂ }  — interpolated from CXRO database
const f2_data = {
  // Tin (Z=50) — dominant EUV absorber; Sn 4d edge near 24-29 eV, strong above 60 eV
  Sn: {
    eV: [80, 82, 84, 86, 88, 90, 91, 91.8, 92, 93, 94, 95, 96, 98, 100, 102, 105],
    f2: [4.42, 4.58, 4.75, 4.91, 5.08, 5.24, 5.32, 5.38, 5.40, 5.47, 5.54, 5.61, 5.68, 5.81, 5.94, 6.07, 6.26],
    Z: 50, A: 118.71
  },
  // Oxygen (Z=8)
  O: {
    eV: [80, 82, 84, 86, 88, 90, 91, 91.8, 92, 93, 94, 95, 96, 98, 100, 102, 105],
    f2: [0.632, 0.643, 0.653, 0.664, 0.675, 0.685, 0.691, 0.695, 0.696, 0.701, 0.706, 0.712, 0.717, 0.727, 0.738, 0.748, 0.764],
    Z: 8, A: 16.00
  },
  // Carbon (Z=6)
  C: {
    eV: [80, 82, 84, 86, 88, 90, 91, 91.8, 92, 93, 94, 95, 96, 98, 100, 102, 105],
    f2: [0.294, 0.299, 0.304, 0.309, 0.314, 0.319, 0.322, 0.324, 0.324, 0.327, 0.329, 0.332, 0.334, 0.339, 0.344, 0.349, 0.357],
    Z: 6, A: 12.01
  },
  // Hydrogen (Z=1)
  H: {
    eV: [80, 82, 84, 86, 88, 90, 91, 91.8, 92, 93, 94, 95, 96, 98, 100, 102, 105],
    f2: [0.00118, 0.00116, 0.00114, 0.00112, 0.00110, 0.00108, 0.00107, 0.00107, 0.00107, 0.00106, 0.00105, 0.00104, 0.00103, 0.00101, 0.000994, 0.000977, 0.000952],
    Z: 1, A: 1.008
  },
  // Nitrogen (Z=7) — for amino/cyano ligands
  N: {
    eV: [80, 82, 84, 86, 88, 90, 91, 91.8, 92, 93, 94, 95, 96, 98, 100, 102, 105],
    f2: [0.451, 0.459, 0.466, 0.474, 0.482, 0.489, 0.493, 0.496, 0.497, 0.501, 0.505, 0.509, 0.512, 0.520, 0.528, 0.536, 0.548],
    Z: 7, A: 14.01
  },
  // Fluorine (Z=9) — for fluorinated ligands
  F: {
    eV: [80, 82, 84, 86, 88, 90, 91, 91.8, 92, 93, 94, 95, 96, 98, 100, 102, 105],
    f2: [0.832, 0.846, 0.860, 0.874, 0.888, 0.902, 0.909, 0.914, 0.916, 0.923, 0.930, 0.937, 0.944, 0.958, 0.972, 0.986, 1.007],
    Z: 9, A: 19.00
  }
};

// Linear interpolation helper
function lerp(xs, ys, x) {
  if (x <= xs[0]) return ys[0];
  if (x >= xs[xs.length - 1]) return ys[ys.length - 1];
  for (let i = 0; i < xs.length - 1; i++) {
    if (x >= xs[i] && x <= xs[i + 1]) {
      const t = (x - xs[i]) / (xs[i + 1] - xs[i]);
      return ys[i] + t * (ys[i + 1] - ys[i]);
    }
  }
  return ys[ys.length - 1];
}

// Get f₂ for element at energy E
function getF2(element, E_eV) {
  const d = f2_data[element];
  if (!d) return 0;
  return lerp(d.eV, d.f2, E_eV);
}

// Atomic photoabsorption cross-section: σ_a = 2·r_e·λ·f₂
function atomicCrossSection(element, E_eV) {
  const lambda_cm = hc / E_eV;
  const f2 = getF2(element, E_eV);
  return 2 * r_e * lambda_cm * f2; // cm²/atom
}

// ════════════════════════════════════════════════════════════════
// MOR PR Molecular Composition Models
// ════════════════════════════════════════════════════════════════

// Predefined MOR architectures
const MOR_PRESETS = {
  "Sn12-butyl": {
    name: "Sn₁₂ Cage (n-Butyl)",
    desc: "[(BuSn)₁₂O₁₄(OH)₆]²⁺ — Inpria-type",
    atoms: (N) => ({ Sn: N, O: Math.round(N * 14/12) + Math.round(N * 6/12), C: N * 4, H: N * 9 + Math.round(N * 6/12) }),
    density: 2.1, // g/cm³ (typical MOR film)
    molRadius: (N) => 0.28 * Math.pow(N, 1/3) + 0.3, // nm
  },
  "Sn12-methyl": {
    name: "Sn₁₂ Cage (Methyl)",
    desc: "[(MeSn)₁₂O₁₄(OH)₆]²⁺",
    atoms: (N) => ({ Sn: N, O: Math.round(N * 14/12) + Math.round(N * 6/12), C: N, H: N * 3 + Math.round(N * 6/12) }),
    density: 2.5,
    molRadius: (N) => 0.25 * Math.pow(N, 1/3) + 0.25,
  },
  "Sn4-cluster": {
    name: "Sn₄ Cluster",
    desc: "Tetranuclear tin-oxo cluster with carboxylate",
    atoms: (N) => ({ Sn: N, O: N * 3, C: N * 6, H: N * 10 }),
    density: 1.8,
    molRadius: (N) => 0.35 * Math.pow(N, 1/3) + 0.3,
  },
  "TTC-macrocycle": {
    name: "TTC Macrocycle",
    desc: "Trinuclear Sn macrocycle (no Sn-OH)",
    atoms: (N) => ({ Sn: N, O: N * 3, C: N * 9, H: N * 8, N: N }),
    density: 1.7,
    molRadius: (N) => 0.4 * Math.pow(N, 1/3) + 0.2,
  },
  "custom": {
    name: "Custom Composition",
    desc: "User-defined atomic composition",
    atoms: () => ({}),
    density: 2.0,
    molRadius: () => 1.0,
  }
};

// Calculate molecular absorption coefficient spectrum
function calcMolecularAbsorption(atoms, density, energies) {
  // μ(E) = (ρ · N_A / M) · Σᵢ nᵢ · σᵢ(E)
  // where M = molecular weight, nᵢ = # of atoms of element i
  let M = 0;
  for (const [el, n] of Object.entries(atoms)) {
    if (f2_data[el]) M += n * f2_data[el].A;
  }
  if (M === 0) return energies.map(E => ({ E, mu: 0, muPerUm: 0 }));

  const prefactor = density * N_A / M; // atoms/cm³ (molecular number density × atoms/molecule combined)

  return energies.map(E => {
    let sigma_total = 0;
    const contributions = {};
    for (const [el, n] of Object.entries(atoms)) {
      if (f2_data[el] && n > 0) {
        const sig = atomicCrossSection(el, E);
        contributions[el] = n * sig * prefactor;
        sigma_total += n * sig;
      }
    }
    const mu = sigma_total * prefactor; // cm⁻¹
    const muPerUm = mu * 1e-4; // μm⁻¹
    return { E: +E.toFixed(1), mu: +mu.toFixed(0), muPerUm: +muPerUm.toFixed(2), lambda: +(hc / E * 1e7).toFixed(2), ...Object.fromEntries(Object.entries(contributions).map(([k,v])=>[`mu_${k}`, +v.toFixed(0)])) };
  });
}

// Calculate absorption with aging effect
function calcAgedAbsorption(baseAtoms, density, energies, f_SnC, f_SnOSn) {
  // As Sn-C bonds break: each broken Sn-C → Sn-OH (lose CₙH₂ₙ₊₁, gain OH)
  // As Sn-O-Sn hydrolyzes: each broken → 2 Sn-OH (gain H₂O)
  const nSn = baseAtoms.Sn || 0;
  const nC_per_ligand = (baseAtoms.C || 0) / Math.max(1, nSn); // C atoms per Sn-C bond
  const nH_per_ligand = Math.max(0, ((baseAtoms.H || 0) - (baseAtoms.O || 0) * 0) / Math.max(1, nSn));

  const nBrokenSnC = Math.round(nSn * (1 - f_SnC));
  const nBrokenSnOSn = Math.round((baseAtoms.O || 0) * 0.7 * (1 - f_SnOSn)); // 70% of O are in Sn-O-Sn bridges

  const agedAtoms = { ...baseAtoms };
  // Ligand loss: remove C and H, add OH
  agedAtoms.C = Math.max(0, (baseAtoms.C || 0) - nBrokenSnC * nC_per_ligand);
  agedAtoms.H = Math.max(0, (baseAtoms.H || 0) - nBrokenSnC * nH_per_ligand + nBrokenSnC + nBrokenSnOSn * 2);
  agedAtoms.O = (baseAtoms.O || 0) + nBrokenSnC + nBrokenSnOSn; // add OH groups

  // Density change: ligand loss increases density (less organic, more inorganic)
  const organicFraction = f_SnC;
  const agedDensity = density * (1 + 0.3 * (1 - organicFraction)); // up to 30% densification

  return {
    spectrum: calcMolecularAbsorption(agedAtoms, agedDensity, energies),
    atoms: agedAtoms,
    density: agedDensity,
    M: Object.entries(agedAtoms).reduce((s, [el, n]) => s + (f2_data[el] ? n * f2_data[el].A : 0), 0)
  };
}

// Gaussian Process Regression (simplified for calibration)
function gpPredict(trainX, trainY, testX, lengthScale = 5.0, noiseVar = 0.01) {
  const n = trainX.length;
  if (n === 0) return testX.map(() => null);

  // RBF kernel
  const kernel = (x1, x2) => Math.exp(-0.5 * ((x1 - x2) / lengthScale) ** 2);

  // K + σ²I
  const K = [];
  for (let i = 0; i < n; i++) {
    K[i] = [];
    for (let j = 0; j < n; j++) {
      K[i][j] = kernel(trainX[i], trainX[j]) + (i === j ? noiseVar : 0);
    }
  }

  // Invert K (simple for small n)
  try {
    const Kinv = math.inv(K);
    const alpha = math.multiply(Kinv, trainY);

    return testX.map(x => {
      let mean = 0;
      for (let i = 0; i < n; i++) {
        mean += alpha[i] * kernel(x, trainX[i]);
      }
      return mean;
    });
  } catch (e) {
    return testX.map(() => null);
  }
}

// ════════════════════════════════════════════════════════════════
// UI
// ════════════════════════════════════════════════════════════════

const P = {
  bg: "#060b18", s1: "#0c1426", s2: "#131f38", bd: "#1e3058",
  ac: "#00d4ff", euv: "#a855f7", sn: "#fbbf24",
  warn: "#f59e0b", err: "#ef4444", ok: "#22c55e",
  t1: "#e2e8f0", t2: "#94a3b8", t3: "#475569",
  snLine: "#ff6b35", oLine: "#4ecdc4", cLine: "#ffe66d", hLine: "#95e1d3",
};

const S = {
  root: { background: P.bg, minHeight: "100vh", color: P.t1, fontFamily: "'IBM Plex Mono', monospace", fontSize: 12 },
  hdr: { background: `linear-gradient(135deg, ${P.s1}, #080e22)`, borderBottom: `1px solid ${P.bd}`, padding: "14px 20px" },
  title: { fontSize: 16, fontWeight: 800, background: `linear-gradient(90deg, ${P.ac}, ${P.euv})`, WebkitBackgroundClip: "text", WebkitTextFillColor: "transparent" },
  sub: { fontSize: 9, color: P.t3, letterSpacing: "2px", textTransform: "uppercase", marginTop: 2 },
  grid: { display: "grid", gridTemplateColumns: "280px 1fr", minHeight: "calc(100vh - 60px)" },
  side: { background: P.s1, borderRight: `1px solid ${P.bd}`, padding: 12, overflowY: "auto", maxHeight: "calc(100vh - 60px)" },
  main: { padding: 12, overflowY: "auto", maxHeight: "calc(100vh - 60px)" },
  card: { background: P.s2, border: `1px solid ${P.bd}`, borderRadius: 5, padding: 12, marginBottom: 10 },
  secT: { fontSize: 9, fontWeight: 700, textTransform: "uppercase", letterSpacing: "2px", color: P.ac, marginBottom: 6 },
  sl: { width: "100%", accentColor: P.ac, margin: "1px 0" },
  mr: { display: "flex", justifyContent: "space-between", padding: "3px 0", borderBottom: `1px solid ${P.bd}15` },
  ttStyle: { background: P.s1, border: `1px solid ${P.bd}`, fontSize: 10, fontFamily: "inherit" },
  tab: (a) => ({ padding: "6px 10px", fontSize: 10, fontWeight: a ? 700 : 500, color: a ? P.ac : P.t2, background: a ? `${P.ac}11` : "transparent", border: "none", borderBottom: a ? `2px solid ${P.ac}` : "2px solid transparent", cursor: "pointer", fontFamily: "inherit" }),
  sel: { width: "100%", background: P.s1, border: `1px solid ${P.bd}`, borderRadius: 3, padding: "5px 7px", color: P.t1, fontSize: 11, fontFamily: "inherit", outline: "none" },
};

function Sl({ label, unit, value, onChange, min, max, step }) {
  return (
    <div style={{ marginBottom: 5 }}>
      <div style={{ display: "flex", justifyContent: "space-between", fontSize: 10, color: P.t2, marginBottom: 1 }}>
        <span>{label}</span><span style={{ color: P.ac, fontWeight: 600 }}>{typeof value === 'number' ? (value < 0.01 ? value.toExponential(2) : value) : value}{unit && ` ${unit}`}</span>
      </div>
      <input type="range" min={min} max={max} step={step} value={value} onChange={e => onChange(+e.target.value)} style={S.sl} />
    </div>
  );
}

function Mv({ l, v, u, c }) {
  return (
    <div style={S.mr}>
      <span style={{ fontSize: 10, color: P.t2 }}>{l}</span>
      <span style={{ fontSize: 11, fontWeight: 600, color: c || P.t1 }}>{typeof v === 'number' ? (Math.abs(v) > 1000 ? v.toFixed(0) : v < 0.01 ? v.toExponential(2) : v.toFixed(3)) : v} <span style={{ fontSize: 9, color: P.t3 }}>{u}</span></span>
    </div>
  );
}

export default function App() {
  const [preset, setPreset] = useState("Sn12-butyl");
  const [nSn, setNSn] = useState(12);
  const [density, setDensity] = useState(2.1);
  const [filmTh, setFilmTh] = useState(20); // nm
  const [tab, setTab] = useState("spectrum");

  // Custom composition (for custom preset)
  const [customSn, setCustomSn] = useState(12);
  const [customO, setCustomO] = useState(20);
  const [customC, setCustomC] = useState(48);
  const [customH, setCustomH] = useState(114);

  // Aging parameters
  const [f_SnC, setF_SnC] = useState(1.0);
  const [f_SnOSn, setF_SnOSn] = useState(1.0);

  // Experimental calibration points
  const [calPoints, setCalPoints] = useState([]);
  const [calInput, setCalInput] = useState({ E: 91.8, mu: "" });

  // Energy range
  const energies = useMemo(() => {
    const es = [];
    for (let E = 80; E <= 105; E += 0.5) es.push(E);
    return es;
  }, []);

  // Get current composition
  const composition = useMemo(() => {
    const p = MOR_PRESETS[preset];
    if (preset === "custom") {
      return { Sn: customSn, O: customO, C: customC, H: customH };
    }
    return p.atoms(nSn);
  }, [preset, nSn, customSn, customO, customC, customH]);

  const currentDensity = useMemo(() => {
    if (preset === "custom") return density;
    return MOR_PRESETS[preset].density;
  }, [preset, density]);

  // Fresh spectrum
  const freshSpec = useMemo(() => calcMolecularAbsorption(composition, currentDensity, energies), [composition, currentDensity, energies]);

  // Aged spectrum
  const agedResult = useMemo(() => calcAgedAbsorption(composition, currentDensity, energies, f_SnC, f_SnOSn), [composition, currentDensity, energies, f_SnC, f_SnOSn]);

  // GP calibration
  const gpCorrected = useMemo(() => {
    if (calPoints.length < 2) return null;
    const trainX = calPoints.map(p => p.E);
    const trainY = calPoints.map(p => {
      const modelMu = freshSpec.find(s => Math.abs(s.E - p.E) < 0.3);
      return modelMu ? p.mu / modelMu.muPerUm : 1;
    });
    const testX = energies;
    const corrections = gpPredict(trainX, trainY, testX, 5.0, 0.001);
    return freshSpec.map((s, i) => ({
      ...s,
      muPerUm_cal: corrections[i] ? +(s.muPerUm * corrections[i]).toFixed(2) : s.muPerUm,
      correction: corrections[i] ? +corrections[i].toFixed(4) : 1,
    }));
  }, [calPoints, freshSpec, energies]);

  // N-sweep data (μ at 13.5nm vs N_Sn)
  const nSweep = useMemo(() => {
    const result = [];
    for (let n = 1; n <= 20; n++) {
      const p = MOR_PRESETS[preset === "custom" ? "Sn12-butyl" : preset];
      const atoms = p.atoms(n);
      const spec = calcMolecularAbsorption(atoms, p.density, [91.8]);
      const M = Object.entries(atoms).reduce((s, [el, num]) => s + (f2_data[el] ? num * f2_data[el].A : 0), 0);
      result.push({
        N: n,
        mu: spec[0]?.muPerUm || 0,
        M: +M.toFixed(0),
        snFrac: +(atoms.Sn * 118.71 / M * 100).toFixed(1),
        OD: +(spec[0]?.muPerUm * filmTh * 0.001 || 0).toFixed(3), // optical density = μ × d
      });
    }
    return result;
  }, [preset, filmTh]);

  // Combined spectrum data for charting
  const specData = useMemo(() => {
    return freshSpec.map((s, i) => ({
      ...s,
      mu_aged: agedResult.spectrum[i]?.muPerUm || s.muPerUm,
      mu_cal: gpCorrected?.[i]?.muPerUm_cal || null,
    }));
  }, [freshSpec, agedResult, gpCorrected]);

  // Key metrics at 13.5 nm (91.8 eV)
  const at135 = useMemo(() => {
    const fresh = freshSpec.find(s => Math.abs(s.E - 91.8) < 0.3) || { muPerUm: 0 };
    const aged = agedResult.spectrum.find(s => Math.abs(s.E - 91.8) < 0.3) || { muPerUm: 0 };
    const OD_fresh = fresh.muPerUm * filmTh * 0.001;
    const OD_aged = aged.muPerUm * filmTh * 0.001;
    const transmission_fresh = Math.exp(-OD_fresh);
    const transmission_aged = Math.exp(-OD_aged);
    const absorption_fresh = 1 - transmission_fresh;
    const absorption_aged = 1 - transmission_aged;
    return { fresh, aged, OD_fresh, OD_aged, transmission_fresh, transmission_aged, absorption_fresh, absorption_aged };
  }, [freshSpec, agedResult, filmTh]);

  // Molecular weight
  const molWeight = Object.entries(composition).reduce((s, [el, n]) => s + (f2_data[el] ? n * f2_data[el].A : 0), 0);

  // Sn weight fraction
  const snWeightFrac = composition.Sn ? (composition.Sn * 118.71 / molWeight) : 0;

  const addCalPoint = () => {
    const mu = parseFloat(calInput.mu);
    if (!isNaN(mu) && mu > 0) {
      setCalPoints([...calPoints, { E: calInput.E, mu }]);
      setCalInput({ ...calInput, mu: "" });
    }
  };

  const tabs = [["spectrum", "μ(λ) Spectrum"], ["nsweep", "N_Sn Dependence"], ["aging", "Aging Effect"], ["calibrate", "ML Calibration"], ["decomp", "Elemental Decomposition"]];

  return (
    <div style={S.root}>
      <div style={S.hdr}>
        <div style={S.title}>EUV Absorption Coefficient Module — MOR PR</div>
        <div style={S.sub}>CXRO/Henke database · 12–15 nm spectral engine · aging-coupled · GP calibration</div>
      </div>

      <div style={S.grid}>
        <div style={S.side}>
          {/* Material Selection */}
          <div style={{ marginBottom: 10 }}>
            <div style={S.secT}><span style={{ color: P.sn }}>◆</span> Material</div>
            <select value={preset} onChange={e => { setPreset(e.target.value); if (e.target.value !== "custom") setDensity(MOR_PRESETS[e.target.value].density); }} style={S.sel}>
              {Object.entries(MOR_PRESETS).map(([k, v]) => <option key={k} value={k}>{v.name}</option>)}
            </select>
            <div style={{ fontSize: 9, color: P.t3, marginTop: 3 }}>{MOR_PRESETS[preset].desc}</div>
          </div>

          {preset !== "custom" ? (
            <Sl label="Sn atoms (N)" value={nSn} onChange={setNSn} min={1} max={20} step={1} />
          ) : (
            <div>
              <Sl label="Sn" value={customSn} onChange={setCustomSn} min={0} max={30} step={1} />
              <Sl label="O" value={customO} onChange={setCustomO} min={0} max={50} step={1} />
              <Sl label="C" value={customC} onChange={setCustomC} min={0} max={100} step={1} />
              <Sl label="H" value={customH} onChange={setCustomH} min={0} max={200} step={1} />
            </div>
          )}

          <Sl label="Film density" unit="g/cm³" value={currentDensity} onChange={setDensity} min={1.0} max={4.0} step={0.1} />
          <Sl label="Film thickness" unit="nm" value={filmTh} onChange={setFilmTh} min={5} max={50} step={1} />

          {/* Composition summary */}
          <div style={S.card}>
            <div style={{ fontSize: 9, color: P.ac, fontWeight: 700, marginBottom: 4 }}>COMPOSITION</div>
            <Mv l="Formula" v={`Sn${composition.Sn}O${composition.O}C${composition.C}H${composition.H}${composition.N ? `N${composition.N}` : ''}`} />
            <Mv l="Mol. weight" v={molWeight} u="g/mol" />
            <Mv l="Sn wt%" v={(snWeightFrac * 100)} u="%" c={P.sn} />
          </div>

          {/* Aging controls */}
          <div style={{ marginBottom: 10 }}>
            <div style={S.secT}><span style={{ color: P.warn }}>◆</span> Aging State</div>
            <Sl label="f(Sn-C) intact" value={f_SnC} onChange={setF_SnC} min={0} max={1} step={0.01} />
            <Sl label="f(Sn-O-Sn) intact" value={f_SnOSn} onChange={setF_SnOSn} min={0} max={1} step={0.01} />
            <div style={{ fontSize: 9, color: P.t3, lineHeight: 1.5 }}>
              f=1.0: fresh film<br />f=0.5: 50% bonds broken (severe aging)
            </div>
          </div>

          {/* Key metrics */}
          <div style={{ ...S.card, background: `${P.ac}08`, border: `1px solid ${P.ac}22` }}>
            <div style={{ fontSize: 9, color: P.ac, fontWeight: 700, marginBottom: 4 }}>@ 13.5 nm (91.8 eV)</div>
            <Mv l="μ (fresh)" v={at135.fresh.muPerUm} u="μm⁻¹" c={P.ac} />
            <Mv l="μ (aged)" v={at135.aged.muPerUm} u="μm⁻¹" c={f_SnC < 0.9 ? P.warn : P.t1} />
            <Mv l="Δμ/μ₀" v={((at135.aged.muPerUm - at135.fresh.muPerUm) / at135.fresh.muPerUm * 100)} u="%" c={Math.abs(at135.aged.muPerUm - at135.fresh.muPerUm) / at135.fresh.muPerUm > 0.05 ? P.warn : P.ok} />
            <Mv l="OD (fresh)" v={at135.OD_fresh} />
            <Mv l="Absorption %" v={(at135.absorption_fresh * 100)} u="%" c={at135.absorption_fresh > 0.3 ? P.ok : P.warn} />
            <Mv l="Abs. % (aged)" v={(at135.absorption_aged * 100)} u="%" />
          </div>
        </div>

        {/* MAIN */}
        <div style={S.main}>
          <div style={{ display: "flex", gap: 2, marginBottom: 12, borderBottom: `1px solid ${P.bd}` }}>
            {tabs.map(([k, l]) => <button key={k} style={S.tab(tab === k)} onClick={() => setTab(k)}>{l}</button>)}
          </div>

          {/* SPECTRUM TAB */}
          {tab === "spectrum" && (
            <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 10 }}>
              <div style={{ ...S.card, gridColumn: "1/-1" }}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>Linear Absorption Coefficient μ(E) — Fresh vs Aged</div>
                <ResponsiveContainer width="100%" height={300}>
                  <ComposedChart data={specData}>
                    <CartesianGrid strokeDasharray="3 3" stroke={P.bd} />
                    <XAxis dataKey="E" tick={{ fill: P.t2, fontSize: 9 }} label={{ value: "Photon Energy (eV)", position: "insideBottom", offset: -5, fill: P.t3, fontSize: 9 }} />
                    <YAxis tick={{ fill: P.t2, fontSize: 9 }} label={{ value: "μ (μm⁻¹)", angle: -90, position: "insideLeft", fill: P.t3, fontSize: 9 }} />
                    <Tooltip contentStyle={S.ttStyle} />
                    <Line type="monotone" dataKey="muPerUm" stroke={P.ac} name="Fresh" strokeWidth={2.5} dot={false} />
                    <Line type="monotone" dataKey="mu_aged" stroke={P.warn} name="Aged" strokeWidth={2} dot={false} strokeDasharray={f_SnC < 0.99 || f_SnOSn < 0.99 ? "none" : "5 5"} />
                    {gpCorrected && <Line type="monotone" dataKey="mu_cal" stroke={P.ok} name="ML-Calibrated" strokeWidth={2} dot={false} />}
                    <ReferenceLine x={91.8} stroke={P.euv} strokeDasharray="4 4" label={{ value: "13.5nm", fill: P.euv, fontSize: 9 }} />
                    {calPoints.map((p, i) => (
                      <ReferenceLine key={i} x={p.E} stroke={P.sn} strokeDasharray="2 2" />
                    ))}
                    <Legend wrapperStyle={{ fontSize: 9 }} />
                  </ComposedChart>
                </ResponsiveContainer>
              </div>

              <div style={S.card}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>Absorption in Film (λ scale)</div>
                <ResponsiveContainer width="100%" height={220}>
                  <AreaChart data={specData.map(s => ({ ...s, absPercent: +((1 - Math.exp(-s.muPerUm * filmTh * 0.001)) * 100).toFixed(1), absPercent_aged: +((1 - Math.exp(-s.mu_aged * filmTh * 0.001)) * 100).toFixed(1) }))}>
                    <CartesianGrid strokeDasharray="3 3" stroke={P.bd} />
                    <XAxis dataKey="lambda" reversed tick={{ fill: P.t2, fontSize: 9 }} label={{ value: "λ (nm)", position: "insideBottom", offset: -5, fill: P.t3, fontSize: 9 }} />
                    <YAxis tick={{ fill: P.t2, fontSize: 9 }} domain={[0, 100]} label={{ value: "Absorption %", angle: -90, position: "insideLeft", fill: P.t3, fontSize: 9 }} />
                    <Tooltip contentStyle={S.ttStyle} />
                    <Area type="monotone" dataKey="absPercent" stroke={P.ac} fill={P.ac} fillOpacity={0.15} name="Fresh" strokeWidth={2} />
                    <Area type="monotone" dataKey="absPercent_aged" stroke={P.warn} fill={P.warn} fillOpacity={0.08} name="Aged" strokeWidth={1.5} />
                    <ReferenceLine x={13.5} stroke={P.euv} strokeDasharray="4 4" />
                    <Legend wrapperStyle={{ fontSize: 9 }} />
                  </AreaChart>
                </ResponsiveContainer>
              </div>

              <div style={S.card}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>Optical Density (OD = μ·d)</div>
                <ResponsiveContainer width="100%" height={220}>
                  <LineChart data={specData.map(s => ({ E: s.E, OD_fresh: +(s.muPerUm * filmTh * 0.001).toFixed(4), OD_aged: +(s.mu_aged * filmTh * 0.001).toFixed(4) }))}>
                    <CartesianGrid strokeDasharray="3 3" stroke={P.bd} />
                    <XAxis dataKey="E" tick={{ fill: P.t2, fontSize: 9 }} />
                    <YAxis tick={{ fill: P.t2, fontSize: 9 }} />
                    <Tooltip contentStyle={S.ttStyle} />
                    <Line type="monotone" dataKey="OD_fresh" stroke={P.ac} name="OD fresh" strokeWidth={2} dot={false} />
                    <Line type="monotone" dataKey="OD_aged" stroke={P.warn} name="OD aged" strokeWidth={2} dot={false} />
                    <ReferenceLine x={91.8} stroke={P.euv} strokeDasharray="4 4" />
                    <Legend wrapperStyle={{ fontSize: 9 }} />
                  </LineChart>
                </ResponsiveContainer>
              </div>
            </div>
          )}

          {/* N_Sn SWEEP TAB */}
          {tab === "nsweep" && (
            <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 10 }}>
              <div style={S.card}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>μ @ 13.5nm vs N_Sn (at constant architecture)</div>
                <ResponsiveContainer width="100%" height={280}>
                  <ComposedChart data={nSweep}>
                    <CartesianGrid strokeDasharray="3 3" stroke={P.bd} />
                    <XAxis dataKey="N" tick={{ fill: P.t2, fontSize: 9 }} label={{ value: "N (Sn atoms)", position: "insideBottom", offset: -5, fill: P.t3, fontSize: 9 }} />
                    <YAxis tick={{ fill: P.t2, fontSize: 9 }} label={{ value: "μ (μm⁻¹)", angle: -90, position: "insideLeft", fill: P.t3, fontSize: 9 }} />
                    <Tooltip contentStyle={S.ttStyle} />
                    <Bar dataKey="mu" fill={P.ac} fillOpacity={0.3} stroke={P.ac} />
                    <Line type="monotone" dataKey="mu" stroke={P.sn} strokeWidth={2} dot={{ fill: P.sn, r: 3 }} name="μ" />
                    <ReferenceLine x={nSn} stroke={P.euv} strokeDasharray="4 4" label={{ value: "current", fill: P.euv, fontSize: 8 }} />
                  </ComposedChart>
                </ResponsiveContainer>
              </div>

              <div style={S.card}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>Sn Weight Fraction & OD vs N_Sn</div>
                <ResponsiveContainer width="100%" height={280}>
                  <ComposedChart data={nSweep}>
                    <CartesianGrid strokeDasharray="3 3" stroke={P.bd} />
                    <XAxis dataKey="N" tick={{ fill: P.t2, fontSize: 9 }} />
                    <YAxis yAxisId="left" tick={{ fill: P.t2, fontSize: 9 }} label={{ value: "Sn wt%", angle: -90, position: "insideLeft", fill: P.t3, fontSize: 9 }} />
                    <YAxis yAxisId="right" orientation="right" tick={{ fill: P.t2, fontSize: 9 }} label={{ value: "OD", angle: 90, position: "insideRight", fill: P.t3, fontSize: 9 }} />
                    <Tooltip contentStyle={S.ttStyle} />
                    <Line yAxisId="left" type="monotone" dataKey="snFrac" stroke={P.sn} name="Sn wt%" strokeWidth={2} dot={{ fill: P.sn, r: 3 }} />
                    <Line yAxisId="right" type="monotone" dataKey="OD" stroke={P.euv} name="OD" strokeWidth={2} dot={{ fill: P.euv, r: 3 }} />
                    <Legend wrapperStyle={{ fontSize: 9 }} />
                  </ComposedChart>
                </ResponsiveContainer>
              </div>

              <div style={{ ...S.card, gridColumn: "1/-1" }}>
                <div style={S.secT}>Design Insight</div>
                <div style={{ fontSize: 10, color: P.t2, lineHeight: 1.8 }}>
                  {nSweep.find(d => d.N === nSn) && (
                    <>
                      N={nSn}: μ = <b style={{ color: P.ac }}>{nSweep.find(d => d.N === nSn).mu.toFixed(1)} μm⁻¹</b>,
                      Sn wt% = <b style={{ color: P.sn }}>{nSweep.find(d => d.N === nSn).snFrac}%</b>,
                      OD = <b>{nSweep.find(d => d.N === nSn).OD}</b> at {filmTh} nm film
                      <br />
                      {nSweep.find(d => d.N === nSn).OD < 0.1 && <span style={{ color: P.err }}>⚠ OD {"<"} 0.1 — insufficient absorption. Increase film thickness or N_Sn.</span>}
                      {nSweep.find(d => d.N === nSn).OD > 0.5 && <span style={{ color: P.warn }}>⚠ OD {">"} 0.5 — strong absorption gradient through film. Consider thinner film.</span>}
                      {nSweep.find(d => d.N === nSn).OD >= 0.1 && nSweep.find(d => d.N === nSn).OD <= 0.5 && <span style={{ color: P.ok }}>✓ OD in optimal range (0.1–0.5) for uniform exposure.</span>}
                      <br />
                      Inpria-class resist: Sn wt% ≈ 30-45%, μ ≈ 15-25 μm⁻¹. CAR comparison: μ ≈ 4-6 μm⁻¹.
                    </>
                  )}
                </div>
              </div>
            </div>
          )}

          {/* AGING TAB */}
          {tab === "aging" && (
            <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 10 }}>
              <div style={{ ...S.card, gridColumn: "1/-1" }}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>Spectral Change Upon Aging: Δμ(E) = μ_aged(E) − μ_fresh(E)</div>
                <ResponsiveContainer width="100%" height={250}>
                  <ComposedChart data={specData.map(s => ({ E: s.E, delta: +(s.mu_aged - s.muPerUm).toFixed(3), deltaPercent: +(((s.mu_aged - s.muPerUm) / s.muPerUm) * 100).toFixed(2) }))}>
                    <CartesianGrid strokeDasharray="3 3" stroke={P.bd} />
                    <XAxis dataKey="E" tick={{ fill: P.t2, fontSize: 9 }} />
                    <YAxis yAxisId="left" tick={{ fill: P.t2, fontSize: 9 }} label={{ value: "Δμ (μm⁻¹)", angle: -90, position: "insideLeft", fill: P.t3, fontSize: 9 }} />
                    <YAxis yAxisId="right" orientation="right" tick={{ fill: P.t2, fontSize: 9 }} label={{ value: "Δμ/μ₀ (%)", angle: 90, position: "insideRight", fill: P.t3, fontSize: 9 }} />
                    <Tooltip contentStyle={S.ttStyle} />
                    <Area yAxisId="left" type="monotone" dataKey="delta" stroke={P.warn} fill={P.warn} fillOpacity={0.2} name="Δμ" strokeWidth={2} />
                    <Line yAxisId="right" type="monotone" dataKey="deltaPercent" stroke={P.err} name="Δμ/μ₀ %" strokeWidth={1.5} dot={false} />
                    <ReferenceLine x={91.8} stroke={P.euv} strokeDasharray="4 4" />
                    <ReferenceLine yAxisId="left" y={0} stroke={P.t3} />
                    <Legend wrapperStyle={{ fontSize: 9 }} />
                  </ComposedChart>
                </ResponsiveContainer>
              </div>

              <div style={S.card}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>Physical Mechanism</div>
                <div style={{ fontSize: 10, color: P.t2, lineHeight: 1.8 }}>
                  <div>f(Sn-C) = <b style={{ color: f_SnC < 0.9 ? P.warn : P.ok }}>{f_SnC.toFixed(2)}</b> → {((1 - f_SnC) * 100).toFixed(0)}% ligands lost</div>
                  <div>f(Sn-O-Sn) = <b style={{ color: f_SnOSn < 0.9 ? P.warn : P.ok }}>{f_SnOSn.toFixed(2)}</b> → {((1 - f_SnOSn) * 100).toFixed(0)}% backbone hydrolyzed</div>
                  <br />
                  <div style={{ color: P.t3 }}>
                    Ligand loss (Sn-C → Sn-OH):<br />
                    → Removes low-Z atoms (C, H) → <b>μ increases</b> (denser Sn-O)<br />
                    → But reduces photoreactive sites → Φ_eff decreases<br /><br />
                    Backbone hydrolysis (Sn-O-Sn → 2 Sn-OH):<br />
                    → Adds H₂O → density decrease → <b>μ may decrease slightly</b><br />
                    → But fragmentation enables dissolution<br /><br />
                    Net: μ typically <b>increases by 5-20%</b> with aging<br />
                    (fewer light atoms, more compact Sn-O network)
                  </div>
                </div>
              </div>

              <div style={S.card}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>Aged Composition</div>
                <Mv l="Sn" v={agedResult.atoms.Sn} c={P.sn} />
                <Mv l="O" v={agedResult.atoms.O} c={P.oLine} />
                <Mv l="C" v={agedResult.atoms.C} c={P.cLine} />
                <Mv l="H" v={agedResult.atoms.H} c={P.hLine} />
                <Mv l="Aged density" v={agedResult.density} u="g/cm³" />
                <Mv l="Aged M" v={agedResult.M} u="g/mol" />
                <Mv l="Aged Sn wt%" v={(agedResult.atoms.Sn * 118.71 / agedResult.M * 100)} u="%" c={P.sn} />
              </div>
            </div>
          )}

          {/* ML CALIBRATION TAB */}
          {tab === "calibrate" && (
            <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 10 }}>
              <div style={S.card}>
                <div style={S.secT}>Add Experimental μ Data Points</div>
                <div style={{ fontSize: 10, color: P.t2, marginBottom: 8 }}>
                  Enter measured μ values at known photon energies. The Gaussian Process will learn a smooth correction function to map the CXRO model to your experimental data.
                </div>
                <div style={{ display: "flex", gap: 6, marginBottom: 8 }}>
                  <div style={{ flex: 1 }}>
                    <div style={{ fontSize: 9, color: P.t3 }}>Energy (eV)</div>
                    <input type="number" value={calInput.E} onChange={e => setCalInput({ ...calInput, E: +e.target.value })} style={{ ...S.sel, width: "100%" }} step={0.1} />
                  </div>
                  <div style={{ flex: 1 }}>
                    <div style={{ fontSize: 9, color: P.t3 }}>μ_exp (μm⁻¹)</div>
                    <input type="number" value={calInput.mu} onChange={e => setCalInput({ ...calInput, mu: e.target.value })} style={{ ...S.sel, width: "100%" }} placeholder="e.g. 18.5" />
                  </div>
                  <button onClick={addCalPoint} style={{ padding: "4px 10px", background: P.ac, color: P.bg, border: "none", borderRadius: 3, fontWeight: 700, cursor: "pointer", alignSelf: "flex-end", fontFamily: "inherit", fontSize: 11 }}>Add</button>
                </div>

                {calPoints.length > 0 && (
                  <div style={{ marginBottom: 8 }}>
                    <div style={{ fontSize: 9, color: P.t3, marginBottom: 3 }}>Calibration points ({calPoints.length})</div>
                    {calPoints.map((p, i) => (
                      <div key={i} style={{ display: "flex", justifyContent: "space-between", fontSize: 10, color: P.t2, padding: "2px 0" }}>
                        <span>{p.E} eV → μ = {p.mu} μm⁻¹</span>
                        <button onClick={() => setCalPoints(calPoints.filter((_, j) => j !== i))} style={{ background: "none", border: "none", color: P.err, cursor: "pointer", fontSize: 10, fontFamily: "inherit" }}>✕</button>
                      </div>
                    ))}
                  </div>
                )}
                <div style={{ fontSize: 9, color: P.t3, lineHeight: 1.5 }}>
                  Need ≥ 2 points for GP regression. Typical lab data sources: ellipsometry (n,k → μ), transmission measurement, or FTIR-derived estimates.
                </div>
              </div>

              <div style={S.card}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>Model vs Experiment</div>
                {gpCorrected ? (
                  <ResponsiveContainer width="100%" height={280}>
                    <ComposedChart data={gpCorrected}>
                      <CartesianGrid strokeDasharray="3 3" stroke={P.bd} />
                      <XAxis dataKey="E" tick={{ fill: P.t2, fontSize: 9 }} />
                      <YAxis tick={{ fill: P.t2, fontSize: 9 }} />
                      <Tooltip contentStyle={S.ttStyle} />
                      <Line type="monotone" dataKey="muPerUm" stroke={P.ac} name="CXRO Model" strokeWidth={2} dot={false} />
                      <Line type="monotone" dataKey="muPerUm_cal" stroke={P.ok} name="GP-Calibrated" strokeWidth={2.5} dot={false} />
                      {calPoints.map((p, i) => (
                        <ReferenceLine key={i} x={p.E} stroke={P.sn} strokeDasharray="2 2" />
                      ))}
                      <Legend wrapperStyle={{ fontSize: 9 }} />
                    </ComposedChart>
                  </ResponsiveContainer>
                ) : (
                  <div style={{ fontSize: 10, color: P.t3, textAlign: "center", padding: "40px 20px" }}>
                    Add ≥ 2 experimental data points to enable GP regression calibration
                  </div>
                )}
              </div>

              {gpCorrected && (
                <div style={{ ...S.card, gridColumn: "1/-1" }}>
                  <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>GP Correction Factor κ(E) = μ_calibrated / μ_CXRO</div>
                  <ResponsiveContainer width="100%" height={180}>
                    <LineChart data={gpCorrected}>
                      <CartesianGrid strokeDasharray="3 3" stroke={P.bd} />
                      <XAxis dataKey="E" tick={{ fill: P.t2, fontSize: 9 }} />
                      <YAxis tick={{ fill: P.t2, fontSize: 9 }} domain={['auto', 'auto']} />
                      <Tooltip contentStyle={S.ttStyle} />
                      <Line type="monotone" dataKey="correction" stroke={P.euv} name="κ(E)" strokeWidth={2} dot={false} />
                      <ReferenceLine y={1} stroke={P.t3} strokeDasharray="4 4" />
                      <ReferenceLine x={91.8} stroke={P.euv} strokeDasharray="4 4" />
                    </LineChart>
                  </ResponsiveContainer>
                  <div style={{ fontSize: 9, color: P.t3 }}>κ = 1.0 means CXRO model is exact. Deviations indicate near-edge effects, chemical state shifts, or density corrections not captured by atomic cross-sections.</div>
                </div>
              )}
            </div>
          )}

          {/* ELEMENTAL DECOMPOSITION TAB */}
          {tab === "decomp" && (
            <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 10 }}>
              <div style={{ ...S.card, gridColumn: "1/-1" }}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>Elemental Contribution to μ(E)</div>
                <ResponsiveContainer width="100%" height={300}>
                  <AreaChart data={freshSpec}>
                    <CartesianGrid strokeDasharray="3 3" stroke={P.bd} />
                    <XAxis dataKey="E" tick={{ fill: P.t2, fontSize: 9 }} />
                    <YAxis tick={{ fill: P.t2, fontSize: 9 }} />
                    <Tooltip contentStyle={S.ttStyle} />
                    {composition.Sn > 0 && <Area type="monotone" dataKey="mu_Sn" stackId="1" stroke={P.snLine} fill={P.snLine} fillOpacity={0.6} name={`Sn (×${composition.Sn})`} />}
                    {composition.O > 0 && <Area type="monotone" dataKey="mu_O" stackId="1" stroke={P.oLine} fill={P.oLine} fillOpacity={0.5} name={`O (×${composition.O})`} />}
                    {composition.C > 0 && <Area type="monotone" dataKey="mu_C" stackId="1" stroke={P.cLine} fill={P.cLine} fillOpacity={0.4} name={`C (×${composition.C})`} />}
                    <ReferenceLine x={91.8} stroke={P.euv} strokeDasharray="4 4" />
                    <Legend wrapperStyle={{ fontSize: 9 }} />
                  </AreaChart>
                </ResponsiveContainer>
              </div>

              <div style={S.card}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>Contribution @ 91.8 eV</div>
                {(() => {
                  const total = at135.fresh.mu;
                  const contribs = [];
                  for (const [el, n] of Object.entries(composition)) {
                    if (f2_data[el] && n > 0) {
                      const sig = atomicCrossSection(el, 91.8);
                      const c = n * sig * currentDensity * N_A / molWeight;
                      contribs.push({ el, n, contrib: c, pct: c / total * 100 * 1e4 });
                    }
                  }
                  contribs.sort((a, b) => b.pct - a.pct);
                  return contribs.map((c, i) => (
                    <div key={i} style={{ display: "flex", alignItems: "center", gap: 6, marginBottom: 4 }}>
                      <div style={{ width: 25, fontSize: 10, fontWeight: 700, color: { Sn: P.snLine, O: P.oLine, C: P.cLine, H: P.hLine, N: "#ff9ff3", F: "#48dbfb" }[c.el] || P.t1 }}>{c.el}</div>
                      <div style={{ flex: 1, height: 6, background: P.bd, borderRadius: 3, overflow: "hidden" }}>
                        <div style={{ width: `${Math.min(100, c.pct)}%`, height: "100%", borderRadius: 3, background: { Sn: P.snLine, O: P.oLine, C: P.cLine, H: P.hLine }[c.el] || P.t3 }} />
                      </div>
                      <div style={{ width: 45, fontSize: 10, textAlign: "right", color: P.t2 }}>{c.pct.toFixed(1)}%</div>
                    </div>
                  ));
                })()}
                <div style={{ marginTop: 8, fontSize: 9, color: P.t3 }}>
                  Sn dominates EUV absorption due to giant 4d cross-section resonance near 92 eV. This is why MOR outperforms CAR by ~5× in photon capture.
                </div>
              </div>

              <div style={S.card}>
                <div style={{ fontSize: 11, fontWeight: 700, marginBottom: 6 }}>Cross-Section Comparison @ 91.8 eV</div>
                {Object.entries(f2_data).map(([el, d]) => (
                  <div key={el} style={{ display: "flex", justifyContent: "space-between", padding: "3px 0", borderBottom: `1px solid ${P.bd}15` }}>
                    <span style={{ fontSize: 10, fontWeight: 700, color: { Sn: P.snLine, O: P.oLine, C: P.cLine, H: P.hLine, N: "#ff9ff3", F: "#48dbfb" }[el] || P.t1 }}>{el} (Z={d.Z})</span>
                    <span style={{ fontSize: 10, color: P.t2 }}>
                      f₂ = {getF2(el, 91.8).toFixed(3)} → σ = {(atomicCrossSection(el, 91.8) * 1e20).toFixed(2)} × 10⁻²⁰ cm²
                    </span>
                  </div>
                ))}
                <div style={{ marginTop: 8, fontSize: 9, color: P.t3 }}>
                  σ_Sn / σ_C ≈ {(atomicCrossSection("Sn", 91.8) / atomicCrossSection("C", 91.8)).toFixed(0)}× — Sn absorbs ~{(atomicCrossSection("Sn", 91.8) / atomicCrossSection("C", 91.8)).toFixed(0)}× more per atom than carbon at 13.5nm
                </div>
              </div>
            </div>
          )}

        </div>
      </div>
    </div>
  );
}
