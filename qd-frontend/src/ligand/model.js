// Ligand-shell thermodynamics of a dot [MX]_m[ML_n]_k0 from the qd pipeline's props/detachment.json
// (ladder down and up) and props/stability.json (free energies of the dot, bulk MX and monomers).
// Everything is evaluated on the pipeline's temperature grid (50-800 K, 25 K steps), in eV.
//
//   ladder state k (units ML_n on the dot), A_k(T) = G_k - G_k0 - (k - k0) G°_MLn:
//     the cumulative removal free energy (k < k0) or addition free energy (k > k0)
//   grand potential   Omega_k - Omega_k0 = A_k - (k - k0) dmu - T (S_k - S_k0),
//                     dmu = mu_MLn - G°_MLn = kT ln(c_free / 1 M)
//   cost of the unit between k - 1 and k (desorption, positive = bound): A_(k-1) - A_k
//   oleylamine:       ML_n + p OAm <=> ML_n(OAm)_p,  Z_A = 1 + sum_p K_p c_OAm^p
//   closed system:    (k0 - <k>) c_dot + C_add = Z_A c_free
//   formation:        m MX + k0 ML_n <=> [MX]_m[ML_n]_k0 with the monomer as supersaturation against bulk MX,
//                     mu_MX = G_bulk + kT ln S; critical S* where the dot as built is in equilibrium

export const KB = 8.617333262e-5;
export const LN10 = Math.LN10;
export const UNITS = { eV: [1, "eV", 2], kcal: [23.060548, "kcal/mol", 1], kJ: [96.4853321, "kJ/mol", 0] };

const lse = (xs) => {
  const m = Math.max(...xs);
  if (!isFinite(m)) return m;
  let s = 0;
  for (const x of xs) s += Math.exp(x - m);
  return m + Math.log(s);
};

function bisect(f, lo, hi, n = 200) {
  let flo = f(lo);
  for (let i = 0; i < n; i++) {
    const mid = 0.5 * (lo + hi), fm = f(mid);
    if ((fm > 0) === (flo > 0)) { lo = mid; flo = fm; } else hi = mid;
    if (hi - lo < 1e-10) break;
  }
  return 0.5 * (lo + hi);
}

// "CdCl2" -> [{t:"CdCl"},{t:"2",sub:true}]
export const formulaParts = (f) =>
  String(f).split(/(\d+)/).filter(Boolean).map((t) => ({ t, sub: /^\d+$/.test(t) }));

export function buildModel(det, stab) {
  const T = det.map?.T || stab?.temperatures;
  const u = stab?.summary?.units || {};
  const steps = det.steps || [], up = det.steps_up || [];
  const unit = det.summary?.unit || `${u.M}${u.X}${u.q > 1 ? u.q : ""}`;
  const M = u.M || unit.match(/^[A-Z][a-z]?/)[0];
  const L = u.X || unit.replace(/^[A-Z][a-z]?/, "").replace(/\d+$/, "");
  const n = u.q || Number(unit.match(/(\d+)$/)?.[1] || 1);
  const k0 = u.m_MXq ?? null, m = u.n_MA ?? null;
  const hasConf = steps.some((s) => s.S_conf_kB !== undefined);
  const S0 = det.configurational?.intact?.S_conf_kB ?? 0;

  // states: k0 - s (removals) and k0 + j (additions), A_k on T, S_k relative to the intact dot (k_B)
  const zero = T.map(() => 0);
  const states = [{ k: 0, A: zero, S: 0, dE: 0 }];          // offsets dk = k - k0 for now
  let acc = zero, accE = 0;
  for (const s of steps) {
    acc = acc.map((v, i) => v + s.dG_eV[i]); accE += s.dE_eV;
    states.push({ k: -s.k, A: acc, S: s.S_conf_kB ?? 0, dE: accE, step: s });
  }
  acc = zero; accE = 0;
  for (const s of up) {
    acc = acc.map((v, i) => v + s.dG_eV[i]); accE += s.dE_eV;
    states.push({ k: -s.k, A: acc, S: s.S_conf_kB ?? 0, dE: accE, step: s });   // s.k < 0 for additions
  }
  states.sort((a, b) => a.k - b.k);
  for (const st of states) { st.dk = st.k; st.k = (k0 ?? 0) + st.dk; }

  const byDk = new Map(states.map((s) => [s.dk, s]));
  // units between consecutive states, indexed by the upper state (dk from dkMin+1 to dkMax)
  const units = [];
  for (const st of states) {
    const lower = byDk.get(st.dk - 1);
    if (!lower) continue;
    const side = st.dk <= 0 ? "removal" : "addition";
    units.push({
      dk: st.dk, k: st.k, side,
      order: side === "removal" ? 1 - st.dk : st.dk,           // removal 1, 2, ... / addition 1, 2, ...
      cost: T.map((_, i) => lower.A[i] - st.A[i]),
      dE: lower.dE - st.dE,
      dS: lower.S - st.S,                                       // S(k-1) - S(k), k_B: entropy gained on desorption
    });
  }

  const ref = stab?.references || {};
  return {
    T, M, L, n, unit, k0, m, coreX: u.A || null, core: u.M && u.A ? `${u.M}${u.A}` : null,
    states, units, hasUp: up.length > 0, hasConf, S0,
    nDown: steps.length, nUp: up.length,
    capDown: steps.length >= (det._inputs?.settings?.max_levels ?? 12),
    capUp: up.length >= (det.summary?.up?.max_levels ?? det._inputs?.settings?.up ?? 12),
    stopped: det.summary?.stopped || "",
    G_MLn: det.reference?.MXq_G_1M_eV || ref.MXq_monomer?.G_1M_eV || null,
    G_dot: stab?.G_dot_1M_eV || null,
    G_bulk: ref.bulk_MA?.G_per_fu_eV || null,
    G_MX: ref.MA_monomer?.G_1M_eV || null,
    bindPerUnit: stab?.binding_dG_per_unit_eV || null,
    decomposition: stab?.decomposition_dG_eV || null,
    siteClasses: (det.site_classes || []).map((c) => ({ dE: c.dE_eV, mult: c.multiplicity || 1, site: c.site || "",
                                                        molecular: c.molecular })),
  };
}

// ---------------- ladder populations against a reservoir of free ML_n ----------------

// conf: true -> ideal S_conf (upper bound), false -> S_conf = 0 (lower bound)
export function populations(model, iT, dmu, conf) {
  const t = model.T[iT], kT = KB * t;
  const om = model.states.map((s) => s.A[iT] - s.dk * dmu - (conf ? kT * s.S : 0));
  const lw = om.map((o) => -o / kT), z = lse(lw);
  const p = lw.map((l) => Math.exp(l - z));
  const kmean = model.states.reduce((a, s, j) => a + s.k * p[j], 0);
  return { omega: om, p, kmean, dkmean: kmean - (model.k0 ?? 0) };
}

export const dmuOf = (cFree, t) => KB * t * Math.log(Math.max(cFree, 1e-300));

// ---------------- oleylamine complexation ----------------

// lnZ_A = ln(1 + sum_p K_p c^p), K_p = exp(-(dG_p + p*delta)/kT); complexes: OAM_COMPLEXES.complexes[unit]
export function amineLnZ(complexes, iT, t, cOAm, delta = 0) {
  if (!complexes?.length || !(cOAm > 0)) return 0;
  const kT = KB * t;
  const terms = [0, ...complexes.map((c) => -(c.dG_eV[iT] + c.p * delta) / kT + c.p * Math.log(cOAm))];
  return lse(terms);
}

export function amineSpeciation(complexes, iT, t, cOAm, delta = 0) {
  const lnZ = amineLnZ(complexes, iT, t, cOAm, delta), kT = KB * t;
  const out = [{ label: "free", frac: Math.exp(-lnZ) }];
  for (const c of complexes || []) {
    const l = cOAm > 0 ? -(c.dG_eV[iT] + c.p * delta) / kT + c.p * Math.log(cOAm) : -Infinity;
    out.push({ label: c.label, p: c.p, frac: Math.exp(l - lnZ) });
  }
  return { lnZ, species: out };
}

// ---------------- closed system: dots as built, nothing but what they release ----------------

// (k0 - <k>) c_dot + C_add = Z_A c_free; returns c_free, <k>, dmu
// kRef: the shell the dots start from (k0 for dots as built, the previous shell after a wash)
export function solveClosed(model, iT, cDot, cAdd, lnZA, conf, kRef = null) {
  const t = model.T[iT], kT = KB * t, k0 = kRef ?? model.k0 ?? 0;
  const kMin = model.states[0].k;
  const f = (x) => {
    const dmu = kT * x;
    const { kmean } = populations(model, iT, dmu, conf);
    const lhs = (k0 - kmean) * cDot + cAdd;
    return Math.exp(x + lnZA) - Math.max(lhs, 0);
  };
  const hi = Math.log(Math.max(cAdd + (k0 - kMin) * cDot, 1e-300)) - lnZA + 1;
  const x = bisect(f, -700, hi);
  const dmu = kT * x;
  return { cFree: Math.exp(x), dmu, ...populations(model, iT, dmu, conf) };
}

// ---------------- formation of the dot from MX monomers and ML_n ----------------

export const hasFormation = (model) => !!(model.G_dot && model.G_bulk && model.G_MLn && model.m && model.k0 !== null);

// log10 of the critical supersaturation of [MX]_m[ML_n]_k0 at free ML_n concentration cFree:
// G_dot + kT ln c_dot - T S_conf - m (G_bulk + kT ln S) - k0 mu_MLn = 0
export function log10Sstar(model, iT, cFree, cDot, conf) {
  const t = model.T[iT], kT = KB * t;
  const mu = model.G_MLn[iT] + kT * Math.log(cFree);
  const g = model.G_dot[iT] + kT * Math.log(cDot) - (conf ? kT * model.S0 : 0)
          - model.m * model.G_bulk[iT] - model.k0 * mu;
  return g / (model.m * kT * LN10);
}
