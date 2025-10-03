#!/usr/bin/env python3
"""
Interactive Plotly visualizer for Fuzzy + PDOS + COOP exports (lightweight).

- Single heatmap (no duplicates); dropdown only restyles zmin.
- Optional --normalize-coop for COOP panel.
- Optional --downsample N to slice the fuzzy grid.
"""

import argparse
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


# ----------------------------- I/O helpers -----------------------------

def load_fuzzy(npz_path):
    d = np.load(npz_path, allow_pickle=True)
    centres = np.asarray(d["centres"], dtype=float)
    Z = np.asarray(d["intensity"], dtype=float)
    tick_positions = np.asarray(d.get("tick_positions", np.arange(Z.shape[1])), dtype=float)
    tick_labels = [str(x) for x in d.get("tick_labels", np.arange(Z.shape[1]))]
    labels = [str(x) for x in d.get("labels", [])]
    extent = np.asarray(
        d.get("extent", [0.0, float(Z.shape[1]-1), float(centres.min()), float(centres.max())]),
        dtype=float
    )
    ewin = np.asarray(d.get("ewin", [float(centres.min()), float(centres.max())]), dtype=float)
    return dict(centres=centres, Z=Z,
                tick_positions=tick_positions, tick_labels=tick_labels,
                labels=labels, extent=extent, ewin=ewin)


def load_pdos_csv(csv_path):
    df = pd.read_csv(csv_path)
    energy = df.iloc[:, 0].to_numpy(dtype=float)    # Energy_eV
    labels = list(df.columns[1:])                   # cumulative columns (Ycum)
    Ycum = df.iloc[:, 1:].to_numpy(dtype=float)     # (nE, nCurves)
    return energy, labels, Ycum


def load_coop_csv(csv_path):
    df = pd.read_csv(csv_path)
    Ener = df.iloc[:, 0].to_numpy(dtype=float)      # MO_Energy_eV
    pairs = list(df.columns[1:])
    values = {p: df[p].to_numpy(dtype=float) for p in pairs}
    return Ener, pairs, values


# ----------------------------- plotting core -----------------------------

def build_combined_figure(
    fuzzy, pdos_energy, pdos_labels, pdos_Ycum,
    coop_energy, coop_pairs, coop_values,
    ef=None, title="Fuzzy Band Map",
    normalize_coop=False, downsample=1
):
    # --- fuzzy inputs ---
    centres = fuzzy["centres"]
    Z = fuzzy["Z"]
    ewin = fuzzy["ewin"]
    kmin, kmax = float(fuzzy["extent"][0]), float(fuzzy["extent"][1])
    nK = Z.shape[1]

    # optional downsample (slice in both axes)
    ds = max(1, int(downsample))
    if ds > 1:
        Z = Z[::ds, ::ds]
        centres = centres[::ds]
        nK = Z.shape[1]

    # k-path axis
    kx = np.linspace(kmin, kmax, nK)

    # tick positions -> path distance if they looked like indices
    tpos = np.asarray(fuzzy["tick_positions"], dtype=float)
    if tpos.size and (tpos.max() <= nK + 1.5) and nK > 1:
        scale = (kmax - kmin) / (nK - 1)
        tpos_plot = kmin + tpos * scale
    else:
        tpos_plot = tpos
    tlabels = fuzzy["tick_labels"]

    # figure scaffold
    fig = make_subplots(
        rows=1, cols=3, shared_yaxes=True,
        column_widths=[0.6, 0.2, 0.2],
        horizontal_spacing=0.06,
        specs=[[{"type": "heatmap"}, {"type": "xy"}, {"type": "xy"}]],
        subplot_titles=(title, "Stacked PDOS", "COOP sticks"),
    )
    fig.update_layout(template="plotly_white",
                      paper_bgcolor="white", plot_bgcolor="white")

    # ---------------- Fuzzy heatmap (single heatmap + dropdown restyle) ----------------
    Zpos = Z[Z > 1e-9]
    vmax = float(np.percentile(Z, 99.9))
    vmin_base = float(np.percentile(Zpos, 5)) if Zpos.size else 1e-6

    Zm = Z.astype(np.float32)
    Zm[Zm <= 0] = np.nan
    Zlog = np.log10(Zm).astype(np.float32)

    # black background behind fuzzy
    fig.add_shape(
        type="rect", xref="x domain", yref="y domain",
        x0=0, x1=1, y0=0, y1=1,
        fillcolor="black", line=dict(width=0), layer="below"
    )

    # single Heatmap trace
    initial_scaled = 1.0
    heat = go.Heatmap(
        z=Zlog, x=kx.astype(np.float32), y=centres.astype(np.float32),
        colorscale="Inferno",
        zmin=float(np.log10(max(vmin_base, vmax / initial_scaled))),
        zmax=float(np.log10(vmax)),
        colorbar=dict(title="log10 Intensity"),
        hovertemplate="k=%{x:.3f} Å⁻¹<br>E=%{y:.3f} eV<br>log10(I)=%{z:.2f}<extra></extra>",
        # small boosts for speed
        hoverongaps=False
    )
    fig.add_trace(heat, row=1, col=1)
    heatmap_idx = len(fig.data) - 1

    # axes + ticks
    fig.update_yaxes(title_text="Energy (eV)", range=[ewin[0], ewin[1]],
                     row=1, col=1, ticks="outside", tickfont=dict(color="black"))
    fig.update_xaxes(title_text="High-Symmetry k-Path", row=1, col=1,
                     ticks="outside", tickfont=dict(color="black"))
    if tpos_plot.size:
        for x in tpos_plot:
            fig.add_vline(x=float(x), line_color="rgba(200,200,200,0.6)",
                          line_width=1, row=1, col=1)
        fig.update_xaxes(
            tickmode="array",
            tickvals=[float(x) for x in tpos_plot],
            ticktext=tlabels,
            row=1, col=1
        )
    if ef is not None:
        fig.add_hline(y=float(ef), line_dash="dash", line_color="white",
                      line_width=2, row=1, col=1)

    # ---------------- PDOS (stacked cumulative) ----------------
    if not np.array_equal(pdos_energy, centres):
        Ycum_use = np.empty((centres.size, pdos_Ycum.shape[1]), dtype=np.float32)
        for j in range(pdos_Ycum.shape[1]):
            Ycum_use[:, j] = np.interp(centres, pdos_energy, pdos_Ycum[:, j]).astype(np.float32)
    else:
        Ycum_use = pdos_Ycum.astype(np.float32)

    palette = (go.Figure().layout.template.layout.colorway
               or ["#636EFA","#EF553B","#00CC96","#AB63FA","#FFA15A",
                   "#19D3F3","#FF6692","#B6E880","#FF97FF","#FECB52"])

    for j, lab in enumerate(pdos_labels):
        ycum = Ycum_use[:, j]
        fig.add_trace(
            go.Scatter(
                x=ycum, y=centres.astype(np.float32),
                mode="lines",
                line=dict(width=0.5, color="rgba(0,0,0,0)"),
                fill="tonextx" if j > 0 else "tozerox",
                fillcolor=palette[j % len(palette)],
                name=lab,
                hovertemplate=f"{lab}: %{x:.3f}<br>E=%{{y:.3f}} eV<extra></extra>",
            ),
            row=1, col=2
        )

    if Ycum_use.shape[1] > 0:
        total = Ycum_use[:, -1]
        fig.add_trace(
            go.Scatter(x=total, y=centres.astype(np.float32), mode="lines",
                       line=dict(color="black", width=2),
                       name="Total DOS",
                       hovertemplate="Total: %{x:.3f}<br>E=%{y:.3f} eV<extra></extra>"),
            row=1, col=2
        )
        xmax = float(max(total.max(), 1e-12)) * 1.05
        fig.update_xaxes(range=[0, xmax], row=1, col=2)

    fig.update_xaxes(title_text="DOS (a.u.)", row=1, col=2)
    fig.update_yaxes(title_text="Energy (eV)", range=[ewin[0], ewin[1]], row=1, col=2)
    if ef is not None:
        fig.add_hline(y=float(ef), line_dash="dash", line_color="black",
                      line_width=1.5, row=1, col=2)

    # ---------------- COOP (interactive sticks; optional normalization) ----------------
    coop_mask = (coop_energy >= ewin[0]) & (coop_energy <= ewin[1])
    Ener = coop_energy[coop_mask]
    if Ener.size == 0 and coop_energy.size:
        e_lo = float(min(ewin[0], np.nanmin(coop_energy)))
        e_hi = float(max(ewin[1], np.nanmax(coop_energy)))
        for c in (1, 2, 3):
            fig.update_yaxes(range=[e_lo, e_hi], row=1, col=c)
        coop_mask = (coop_energy >= e_lo) & (coop_energy <= e_hi)
        Ener = coop_energy[coop_mask]

    if normalize_coop:
        gmax = max(
            (np.nanmax(np.abs(coop_values[p][coop_mask])) for p in coop_pairs if coop_values[p][coop_mask].size),
            default=0.0
        )
        scale = (1.0 / gmax) if gmax > 0 else 1.0
        x_title = "COOP (normalized)"
        x_range = [-1.05, 1.05]
    else:
        scale = 1.0
        vals = []
        for p in coop_pairs:
            a = coop_values[p][coop_mask]
            if a.size:
                vals.append(np.nanmax(np.abs(a)))
        vmax_abs = max(vals) if vals else 1.0
        pad = 0.05 * vmax_abs
        x_title = "COOP (a.u.)"
        x_range = [-(vmax_abs + pad), (vmax_abs + pad)]

    pair_colors = {p: palette[i % len(palette)] for i, p in enumerate(coop_pairs)}

    total_points = sum(int(coop_values[p][coop_mask].size) for p in coop_pairs)
    add_markers = total_points <= 20000

    for p in coop_pairs:
        v = (coop_values[p][coop_mask] * scale).astype(np.float32)
        if v.size == 0:
            continue
        xs, ys = [], []
        for yi, xv in zip(Ener, v):
            xs.extend([0.0, float(xv), None])
            ys.extend([float(yi),  float(yi), None])
        fig.add_trace(
            go.Scattergl(
                x=xs, y=ys, mode="lines",
                line=dict(color=pair_colors[p], width=2),
                name=p,
                hoverinfo="skip",
                showlegend=False
            ),
            row=1, col=3
        )
        if add_markers:
            fig.add_trace(
                go.Scattergl(
                    x=v, y=Ener,
                    mode="markers",
                    marker=dict(size=5, color=pair_colors[p]),
                    name=p,
                    hovertemplate=f"{p}<br>E=%{{y:.3f}} eV<br>COOP=%{{x:.3f}}<extra></extra>",
                ),
                row=1, col=3
            )

    fig.update_xaxes(range=x_range, title_text=x_title, row=1, col=3)
    fig.update_yaxes(title_text="Energy (eV)", row=1, col=3)
    if ef is not None:
        fig.add_hline(y=float(ef), line_dash="dash", line_color="black",
                      line_width=1.5, row=1, col=3)

    # ---------------- Dropdown (updates only zmin on the fuzzy heatmap) ----------------
    scaled_opts = [1, 10, 100, 1000, 10000]
    buttons = []
    for sval in scaled_opts:
        vmin_opt_log10 = float(np.log10(max(vmin_base, vmax / float(sval))))
        buttons.append(dict(
            label=f"scaled_vmin={sval}",
            method="restyle",
            args=[{"zmin": [vmin_opt_log10]}, [heatmap_idx]]
        ))

    fig.update_layout(
        updatemenus=[dict(
            type="dropdown",
            buttons=buttons,
            direction="down",
            showactive=True,
            x=0.98, y=1.10, xanchor="right", yanchor="top",
            pad={"r": 4, "t": 4}
        )],
        legend=dict(
            orientation="h",
            x=0.0,  xanchor="left",
            y=1.10, yanchor="top"
        ),
        margin=dict(l=80, r=40, t=92, b=70),
        height=800,
    )

    return fig


# ----------------------------- CLI -----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fuzzy", required=True, help="Path to fuzzy_data.npz")
    ap.add_argument("--pdos",  required=True, help="Path to pdos_data.csv")
    ap.add_argument("--coop",  required=True, help="Path to coop_data.csv")
    ap.add_argument("--out",   default="fuzzy_pdos_coop_interactive.html", help="Output HTML")
    ap.add_argument("--normalize-coop", action="store_true", help="Normalize COOP to [-1,1]")
    ap.add_argument("--ef", type=float, default=None, help="Fermi/midgap energy for dashed line")
    ap.add_argument("--title", type=str, default="Fuzzy Band Map", help="Title for fuzzy panel")
    ap.add_argument("--downsample", type=int, default=1, help="Slice the fuzzy grid by this factor (1=no downsample)")
    args = ap.parse_args()

    fuzzy = load_fuzzy(args.fuzzy)
    pdos_energy, pdos_labels, pdos_Ycum = load_pdos_csv(args.pdos)
    coop_energy, coop_pairs, coop_values = load_coop_csv(args.coop)

    fig = build_combined_figure(
        fuzzy, pdos_energy, pdos_labels, pdos_Ycum,
        coop_energy, coop_pairs, coop_values,
        ef=args.ef, title=args.title,
        normalize_coop=bool(args.normalize_coop),
        downsample=args.downsample
    )
    fig.write_html(args.out, include_plotlyjs="cdn", full_html=True)
    print(f"✓ Wrote interactive HTML → {args.out}")

if __name__ == "__main__":
    main()

