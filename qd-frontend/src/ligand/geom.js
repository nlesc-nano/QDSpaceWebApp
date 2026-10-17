// Geometry analysis in the browser, as the pipeline's structure step does it: bonds are opposite-element
// pairs closer than the cutoffs of props/structure.json (cutoffs_A, e.g. {"Cd-Se": 3.19, "Cd-Cl": 3.08}).

export function parseXyz(text) {
  const lines = String(text).replace(/\r\n/g, "\n").split("\n");
  const n = parseInt(lines[0], 10);
  const sym = [], pos = [];
  for (let i = 2; i < 2 + n; i++) {
    const t = lines[i].trim().split(/\s+/);
    sym.push(t[0]);
    pos.push([+t[1], +t[2], +t[3]]);
  }
  return { sym, pos };
}

const dist = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);

// bonds by pair: {"Cd-Se": [{d, core}], ...}; core: both atoms have role "core" (roles in the same order)
export function bonds(g, cutoffs, roles = null) {
  const pairs = Object.entries(cutoffs).map(([p, c]) => [p, p.split("-"), c]);
  const out = Object.fromEntries(pairs.map(([p]) => [p, []]));
  const nb = g.sym.map(() => ({}));
  for (let i = 0; i < g.sym.length; i++) {
    for (let j = i + 1; j < g.sym.length; j++) {
      for (const [p, [a, b], c] of pairs) {
        if (!((g.sym[i] === a && g.sym[j] === b) || (g.sym[i] === b && g.sym[j] === a))) continue;
        const d = dist(g.pos[i], g.pos[j]);
        if (d < c) {
          out[p].push({ d, core: roles ? roles[i] === "core" && roles[j] === "core" : null });
          nb[i][g.sym[j]] = (nb[i][g.sym[j]] || 0) + 1;
          nb[j][g.sym[i]] = (nb[j][g.sym[i]] || 0) + 1;
        }
      }
    }
  }
  return { bonds: out, neighbours: nb };
}

// coordination numbers: {element: {cn: count}}
export function cnCounts(g, neighbours) {
  const out = {};
  g.sym.forEach((e, i) => {
    const cn = Object.values(neighbours[i]).reduce((a, b) => a + b, 0);
    out[e] = out[e] || {};
    out[e][cn] = (out[e][cn] || 0) + 1;
  });
  return out;
}

// per-atom displacement between two geometries of the same atoms after optimal superposition (Kabsch)
export function displacements(a, b) {
  const n = a.pos.length;
  const ca = [0, 1, 2].map((k) => a.pos.reduce((s, p) => s + p[k], 0) / n);
  const cb = [0, 1, 2].map((k) => b.pos.reduce((s, p) => s + p[k], 0) / n);
  const P = a.pos.map((p) => p.map((v, k) => v - ca[k])), Q = b.pos.map((p) => p.map((v, k) => v - cb[k]));
  // H = P^T Q; rotation from the polar decomposition, by quaternion (Horn 1987)
  const H = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < n; i++) for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) H[r][c] += P[i][r] * Q[i][c];
  const [[sxx, sxy, sxz], [syx, syy, syz], [szx, szy, szz]] = H;
  const N = [[sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
             [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
             [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
             [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz]];
  let q = [1, 0, 0, 0];                                     // largest eigenvector by power iteration
  const shift = Math.abs(N.flat().reduce((m, v) => Math.max(m, Math.abs(v)), 0)) * 4;
  for (let it = 0; it < 500; it++) {
    const v = N.map((row, r) => row.reduce((s, x, c) => s + x * q[c], 0) + shift * q[r]);
    const nv = Math.hypot(...v);
    q = v.map((x) => x / nv);
  }
  const [w, x, y, z] = q;
  const R = [[w * w + x * x - y * y - z * z, 2 * (x * y - w * z), 2 * (x * z + w * y)],
             [2 * (x * y + w * z), w * w - x * x + y * y - z * z, 2 * (y * z - w * x)],
             [2 * (x * z - w * y), 2 * (y * z + w * x), w * w - x * x - y * y + z * z]];
  return P.map((p, i) => {
    const rp = [0, 1, 2].map((r) => R[r][0] * p[0] + R[r][1] * p[1] + R[r][2] * p[2]);
    return dist(rp, Q[i]);
  });
}

export const sameAtoms = (a, b) => a && b && a.sym.length === b.sym.length && a.sym.every((s, i) => s === b.sym[i]);

export function histogram(values, lo, hi, width) {
  const n = Math.round((hi - lo) / width), counts = new Array(n).fill(0);
  for (const v of values) {
    const i = Math.floor((v - lo) / width);
    if (i >= 0 && i < n) counts[i]++;
  }
  return counts.map((c, i) => ({ x: lo + (i + 0.5) * width, c }));
}
