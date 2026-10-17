// Re-broadening of a curve stored with Gaussian width sigma0 on a uniform grid: a Gaussian of width s
// convolved with one of width sigma0 is a Gaussian of width sqrt(sigma0² + s²), so any sigma >= sigma0
// follows exactly from the stored curve, with no need for the individual peaks.
//   x0, dx: the grid (x_i = x0 + i dx); y: the stored curve; returns [[x, y]] at width sigma, padded by
//   4 sigma at both ends so the tails are not cut, and with the area kept.
export function rebroaden(x0, dx, y, sigma0, sigma) {
  const s = Math.sqrt(Math.max(0, sigma * sigma - sigma0 * sigma0));
  const n = Math.ceil((4 * s) / dx);
  if (n < 1) return y.map((v, i) => [x0 + i * dx, v]);
  const k = [];
  for (let j = -n; j <= n; j++) k.push(Math.exp(-0.5 * ((j * dx) / s) ** 2));
  const ks = k.reduce((a, b) => a + b, 0);
  const out = [];
  for (let i = -n; i < y.length + n; i++) {
    let v = 0;
    for (let j = -n; j <= n; j++) {
      const m = i - j;
      if (m >= 0 && m < y.length) v += y[m] * k[j + n];
    }
    out.push([x0 + i * dx, v / ks]);
  }
  return out;
}
