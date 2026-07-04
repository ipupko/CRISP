#!/usr/bin/env python3
# ##############################################################
#   .oooooo.   ooooooooo.   ooooo  .oooooo..o ooooooooo.
#  d8P'  `Y8b  `888   `Y88. `888' d8P'    `Y8 `888   `Y88.
# 888           888   .d88'  888  Y88bo.       888   .d88'
# 888           888ooo88P'   888   `"Y8888o.   888ooo88P'
# 888           888`88b.     888       `"Y88b  888
# `88b    ooo   888  `88b.   888  oo     .d8P  888
#  `Y8bood8P'  o888o  o888o o888o 8""88888P'  o888o
# ##############################################################
#
# Script : scripts/plot_pca_pass1_3d.py
# Part of: Compass Genomics suite — github.com/ipupko/CRISP
# Version: v0.5.0
#
# Purpose:
#   Interactive 3D PCA scatter — CRISP Step 6 PCA Pass 1.
#   Standalone Plotly HTML with JS-powered controls.
#
#   v1.4.0 changes:
#     - Tooltip annotations on every button (title attribute)
#     - Isolate works for multi-selection (not just single group)
#     - Sub-pop colouring logic:
#         single country selected → sub-pops get distinct COLOURS
#         multiple countries selected → sub-pops use SHAPES within
#                                       each country's colour
#
##########################################################################

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import plotly.graph_objects as go
except ImportError:
    print("[CRISP] ERROR: plotly not found. pip install plotly", file=sys.stderr)
    sys.exit(1)


##########################################################################
# ARGUMENT PARSING
##########################################################################

def parse_args():
    p = argparse.ArgumentParser(
        description="CRISP Step 6 — interactive 3D PCA plot v1.4"
    )
    p.add_argument("--eigenvec",   required=True)
    p.add_argument("--satellite",  required=True)
    p.add_argument("--out-dir",    required=True)
    p.add_argument("--project",    required=True)
    p.add_argument("--pc-x",       type=int, default=1)
    p.add_argument("--pc-y",       type=int, default=2)
    p.add_argument("--pc-z",       type=int, default=3)
    p.add_argument("--admixed",    default="NO")
    p.add_argument("--ellipse-sd", type=float, default=3.0)
    p.add_argument("--static",     default="NO")
    return p.parse_args()


##########################################################################
# PALETTES
##########################################################################

def _pal(var, fallback):
    val = os.environ.get(var, "").strip()
    return val if val else fallback


PALETTES = {
    "STANDARD / LIGHT": {
        "mode":"STANDARD","background":"LIGHT",
        "group_cols":["#1D9E75","#56B4E9","#CC79A7","#FEBC2E",
                      "#8BE0CB","#0072B2","#E76F51","#A8DADC",
                      "#F4A261","#264653","#E9C46A","#2A9D8F"],
        "fail":"#FF5F57","subtext":"#8E8E93",
        "text":"#1C1C1E","bg":"#FFFFFF",
        "panel":"#F5F5F7","grid":"#D1D1D6",
    },
    "COLOURBLIND / LIGHT": {
        "mode":"COLOURBLIND","background":"LIGHT",
        "group_cols":["#0072B2","#56B4E9","#CC79A7","#E69F00",
                      "#F0E442","#009E73","#D55E00","#999999",
                      "#44AA99","#882255","#DDCC77","#117733"],
        "fail":"#D55E00","subtext":"#8E8E93",
        "text":"#1C1C1E","bg":"#FFFFFF",
        "panel":"#F5F5F7","grid":"#D1D1D6",
    },
    "NIGHT / DARK": {
        "mode":"NIGHT","background":"DARK",
        "group_cols":["#2ED9A3","#7FD1FF","#D8B4FE","#FFC857",
                      "#A8F0D8","#58A6FF","#FF9F7F","#B5EAD7",
                      "#FFDAC1","#C7CEEA","#E2F0CB","#FFD1DC"],
        "fail":"#FF7B72","subtext":"#8B949E",
        "text":"#E6EDF3","bg":"#0D1117",
        "panel":"#161B22","grid":"#30363D",
    },
}


def active_palette():
    mode = os.environ.get("CRISP_PAL_MODE","STANDARD")
    bg   = os.environ.get("CRISP_PAL_BACKGROUND","LIGHT")
    return PALETTES.get(f"{mode} / {bg}", PALETTES["STANDARD / LIGHT"])


##########################################################################
# DATA LOADING
##########################################################################

def read_eigenvec(path):
    with open(path) as f:
        first = f.readline().strip()
    if first.startswith("#"):
        df = pd.read_csv(path, sep=r"\s+")
        df.columns = [c.lstrip("#") for c in df.columns]
    else:
        df = pd.read_csv(path, sep=r"\s+", header=None)
        n = df.shape[1]-2
        df.columns = ["FID","IID"]+[f"PC{i+1}" for i in range(n)]
    return df


def load_data(args):
    eigenvec  = read_eigenvec(args.eigenvec)
    satellite = pd.read_csv(args.satellite, sep="\t")
    pc_cols   = [f"PC{i}" for i in sorted({args.pc_x,args.pc_y,args.pc_z})]
    missing   = [c for c in pc_cols if c not in eigenvec.columns]
    if missing:
        print(f"[CRISP] ERROR: PCs not found: {missing}", file=sys.stderr)
        sys.exit(1)

    sat_cols = ["FID","IID","SUPERPOP","OUTLIER","MEMBERSHIP_P"]
    if "SUBPOP" in satellite.columns:
        sat_cols.append("SUBPOP")

    df = satellite[sat_cols].merge(
        eigenvec[["IID"]+pc_cols], on="IID", how="inner"
    )
    df["SUPERPOP"]     = df["SUPERPOP"].astype(str)
    df["OUTLIER"]      = df["OUTLIER"].astype(str)
    df["MEMBERSHIP_P"] = pd.to_numeric(df["MEMBERSHIP_P"], errors="coerce")

    if "SUBPOP" in df.columns:
        df["SUBPOP"] = df["SUBPOP"].astype(str).replace(
            {"—":"","nan":"","None":"","NA":""}
        )
        df.loc[df["SUBPOP"].str.len()==0,"SUBPOP"] = ""
    else:
        df["SUBPOP"] = ""

    for c in pc_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=pc_cols)

    print(f"[CRISP] Samples: {len(df)}")
    has_sub = (df["SUBPOP"].str.len()>0).sum()
    print(f"[CRISP] Samples with SUBPOP: {has_sub}")
    return df


##########################################################################
# ELLIPSOID
##########################################################################

def ellipsoid_surface(centroid, cov, sd, n=25):
    try:
        vals,vecs = np.linalg.eigh(cov)
        vals = np.clip(vals,0,None)
        order = vals.argsort()[::-1]
        vals,vecs = vals[order],vecs[:,order]
        u = np.linspace(0,2*np.pi,n)
        v = np.linspace(0,np.pi,n)
        sx = np.outer(np.cos(u),np.sin(v))
        sy = np.outer(np.sin(u),np.sin(v))
        sz = np.outer(np.ones(n),np.cos(v))
        radii = sd*np.sqrt(vals)
        E = np.zeros((3,n,n))
        for i in range(3):
            E += radii[i]*vecs[:,i:i+1,np.newaxis]*np.array([sx,sy,sz])[i]
        return E[0]+centroid[0],E[1]+centroid[1],E[2]+centroid[2]
    except Exception:
        return None,None,None


##########################################################################
# BUILD PLOTLY FIGURE
##########################################################################

def build_figure(df, args):
    pc_x,pc_y,pc_z = f"PC{args.pc_x}",f"PC{args.pc_y}",f"PC{args.pc_z}"
    pal       = active_palette()
    seq       = pal["group_cols"]

    all_groups = sorted(df["SUPERPOP"].unique().tolist())
    n_grp      = len(all_groups)
    n_out      = int((df["OUTLIER"]=="YES").sum())
    n_tot      = len(df)

    grp_cols = {g: seq[i%len(seq)] for i,g in enumerate(all_groups)}

    # Sub-population map
    subpop_map = {}
    for grp in all_groups:
        subs = sorted(df[(df["SUPERPOP"]==grp)&
                         (df["SUBPOP"].str.len()>0)]["SUBPOP"].unique().tolist())
        if subs:
            subpop_map[grp] = subs

    has_subpops = len(subpop_map) > 0
    if has_subpops:
        for grp,subs in subpop_map.items():
            print(f"[CRISP]   {grp}: {subs}")

    fig = go.Figure()
    trace_meta = []

    # Shapes for multi-country sub-pop mode
    # WHY: When multiple countries are selected, colour = country,
    # shape = sub-population within that country. When only one country
    # is selected, colour = sub-population (much easier to distinguish).
    SUBPOP_SHAPES = ["circle","diamond","square","cross","x",
                     "triangle-up","triangle-down","pentagon","hexagon","star"]

    # ── Group scatter traces ──────────────────────────────────────────
    for grp in all_groups:
        gd   = df[(df["SUPERPOP"]==grp)&(df["OUTLIER"]!="YES")]
        gcol = grp_cols[grp]
        if len(gd)==0: continue

        hover = [
            f"<b>{r['IID']}</b><br>"
            f"Group: {r['SUPERPOP']}<br>"
            f"Sub-pop: {r['SUBPOP'] if r['SUBPOP'] else '—'}<br>"
            f"PC{args.pc_x}: {r[pc_x]:.6f}<br>"
            f"PC{args.pc_y}: {r[pc_y]:.6f}<br>"
            f"PC{args.pc_z}: {r[pc_z]:.6f}<br>"
            f"Membership P: {r['MEMBERSHIP_P']:.4f}"
            for _,r in gd.iterrows()
        ]
        fig.add_trace(go.Scatter3d(
            x=gd[pc_x],y=gd[pc_y],z=gd[pc_z],
            mode="markers",
            name=f"{grp} (n={len(gd)})",
            marker=dict(size=3,color=gcol,opacity=0.70,line=dict(width=0)),
            text=hover,hovertemplate="%{text}<extra></extra>",
            legendgroup=grp,
        ))
        trace_meta.append({"type":"scatter","group":grp,"subpop":"",
                            "base_color":gcol})

        # Ellipsoid
        if args.admixed.upper()=="YES" and len(gd)>=5:
            X3  = gd[[pc_x,pc_y,pc_z]].values
            ex,ey,ez = ellipsoid_surface(X3.mean(axis=0),
                                          np.cov(X3,rowvar=False),
                                          args.ellipse_sd)
            if ex is not None:
                fig.add_trace(go.Surface(
                    x=ex,y=ey,z=ez,opacity=0.08,
                    colorscale=[[0,gcol],[1,gcol]],
                    showscale=False,showlegend=False,
                    hoverinfo="skip",name=f"{grp}_ell",
                    contours=dict(
                        x=dict(highlight=False),
                        y=dict(highlight=False),
                        z=dict(highlight=False),
                    ),
                ))
                trace_meta.append({"type":"ellipsoid","group":grp,"subpop":"",
                                    "base_color":gcol})

    # ── Sub-population traces ─────────────────────────────────────────
    # Two sets of sub-pop traces:
    #   _colour  — used when exactly ONE country selected (colour mode)
    #   _shape   — used when MULTIPLE countries selected (shape mode)
    # WHY two sets: dynamically changing both marker.color and marker.symbol
    # via Plotly.restyle() simultaneously is reliable, but storing both
    # variants as separate traces and toggling visibility is cleaner
    # and avoids race conditions in the JS restyle calls.

    for grp, subs in subpop_map.items():
        gcol = grp_cols[grp]

        for si, sub in enumerate(subs):
            sd_df = df[(df["SUPERPOP"]==grp)&(df["SUBPOP"]==sub)&
                       (df["OUTLIER"]!="YES")]
            if len(sd_df)==0: continue

            # Sub-pop colour (used when single country focused)
            sub_col = seq[(n_grp + si) % len(seq)]

            hover = [
                f"<b>{r['IID']}</b><br>"
                f"Group: {r['SUPERPOP']}<br>"
                f"Sub-pop: <b>{sub}</b><br>"
                f"PC{args.pc_x}: {r[pc_x]:.6f}<br>"
                f"PC{args.pc_y}: {r[pc_y]:.6f}<br>"
                f"PC{args.pc_z}: {r[pc_z]:.6f}"
                for _,r in sd_df.iterrows()
            ]

            # COLOUR mode trace (single-country focus)
            fig.add_trace(go.Scatter3d(
                x=sd_df[pc_x],y=sd_df[pc_y],z=sd_df[pc_z],
                mode="markers",
                name=f"{grp} › {sub} (n={len(sd_df)})",
                marker=dict(
                    size=5, color=sub_col,
                    symbol="circle",
                    opacity=0.88,
                    line=dict(width=1,color=pal["text"]),
                ),
                text=hover,hovertemplate="%{text}<extra></extra>",
                visible=False,
                legendgroup=f"{grp}_sub_col",
                showlegend=False,
            ))
            trace_meta.append({"type":"subpop_col","group":grp,"subpop":sub,
                                "base_color":sub_col,"sub_idx":si})

            # SHAPE mode trace (multi-country focus)
            fig.add_trace(go.Scatter3d(
                x=sd_df[pc_x],y=sd_df[pc_y],z=sd_df[pc_z],
                mode="markers",
                name=f"{grp} › {sub} ◆",
                marker=dict(
                    size=5, color=gcol,
                    symbol=SUBPOP_SHAPES[si%len(SUBPOP_SHAPES)],
                    opacity=0.88,
                    line=dict(width=1,color=pal["text"]),
                ),
                text=hover,hovertemplate="%{text}<extra></extra>",
                visible=False,
                legendgroup=f"{grp}_sub_shp",
                showlegend=False,
            ))
            trace_meta.append({"type":"subpop_shp","group":grp,"subpop":sub,
                                "base_color":gcol,"sub_idx":si})

            # Ellipsoid — single colour set (colour mode only)
            if len(sd_df)>=5:
                X3s = sd_df[[pc_x,pc_y,pc_z]].values
                ex,ey,ez = ellipsoid_surface(
                    X3s.mean(axis=0),
                    np.cov(X3s,rowvar=False),
                    max(1.5,args.ellipse_sd*0.7)
                )
                if ex is not None:
                    fig.add_trace(go.Surface(
                        x=ex,y=ey,z=ez,opacity=0.12,
                        colorscale=[[0,sub_col],[1,sub_col]],
                        showscale=False,showlegend=False,
                        hoverinfo="skip",name=f"{grp}_{sub}_ell",
                        visible=False,
                        contours=dict(
                            x=dict(highlight=False),
                            y=dict(highlight=False),
                            z=dict(highlight=False),
                        ),
                    ))
                    trace_meta.append({"type":"subpop_ell","group":grp,
                                       "subpop":sub,"base_color":sub_col,
                                       "sub_idx":si})

    # ── Outliers ─────────────────────────────────────────────────────
    out = df[df["OUTLIER"]=="YES"]
    if len(out)>0:
        hover_out = [
            f"<b>{r['IID']}</b> ⚠ OUTLIER<br>Group: {r['SUPERPOP']}<br>"
            f"PC{args.pc_x}: {r[pc_x]:.6f}<br>"
            f"PC{args.pc_y}: {r[pc_y]:.6f}<br>"
            f"PC{args.pc_z}: {r[pc_z]:.6f}"
            for _,r in out.iterrows()
        ]
        fig.add_trace(go.Scatter3d(
            x=out[pc_x],y=out[pc_y],z=out[pc_z],
            mode="markers",name=f"Outlier (n={len(out)})",
            marker=dict(size=7,color=pal["fail"],symbol="cross",
                        opacity=0.9,line=dict(width=1,color=pal["fail"])),
            text=hover_out,hovertemplate="%{text}<extra></extra>",
            legendgroup="outlier",
        ))
        trace_meta.append({"type":"outlier","group":"outlier","subpop":"",
                            "base_color":pal["fail"]})

    # ── Centroid diamonds ─────────────────────────────────────────────
    centroids = {}
    for grp in all_groups:
        gd = df[(df["SUPERPOP"]==grp)&(df["OUTLIER"]!="YES")]
        if len(gd)<3: continue
        cx,cy,cz = float(gd[pc_x].mean()),float(gd[pc_y].mean()),float(gd[pc_z].mean())
        centroids[grp] = (cx,cy,cz)
        gcol = grp_cols[grp]
        fig.add_trace(go.Scatter3d(
            x=[cx],y=[cy],z=[cz],
            mode="markers+text",name=f"{grp} ▲",
            text=[grp],textposition="top center",
            textfont=dict(size=9,color=pal["text"]),
            marker=dict(size=7,color=gcol,symbol="diamond",
                        line=dict(width=2,color=pal["text"])),
            showlegend=False,
            hovertemplate=(
                f"<b>{grp} centroid</b><br>"
                f"PC{args.pc_x}: {cx:.6f}<br>"
                f"PC{args.pc_y}: {cy:.6f}<br>"
                f"PC{args.pc_z}: {cz:.6f}<extra></extra>"
            ),
        ))
        trace_meta.append({"type":"centroid","group":grp,"subpop":"",
                            "base_color":gcol})

    n_traces = len(fig.data)

    # ── Palette buttons ───────────────────────────────────────────────
    palette_buttons = []
    for pal_name, p in PALETTES.items():
        seq2   = p["group_cols"]
        g_cols = [seq2[i%len(seq2)] for i in range(n_grp)]
        safe_c = [""]*n_traces
        cscale = [""]*n_traces

        for ti,tm in enumerate(trace_meta):
            if tm["type"] in ("scatter","centroid"):
                gi = all_groups.index(tm["group"]) if tm["group"] in all_groups else 0
                safe_c[ti] = g_cols[gi%len(g_cols)]
            elif tm["type"]=="subpop_shp":
                gi = all_groups.index(tm["group"]) if tm["group"] in all_groups else 0
                safe_c[ti] = g_cols[gi%len(g_cols)]
            elif tm["type"]=="subpop_col":
                safe_c[ti] = seq2[(n_grp + tm.get("sub_idx",0)) % len(seq2)]
            elif tm["type"]=="outlier":
                safe_c[ti] = p["fail"]
            if tm["type"] in ("ellipsoid","subpop_ell"):
                gi = all_groups.index(tm["group"]) if tm["group"] in all_groups else 0
                col = seq2[(n_grp+tm.get("sub_idx",0))%len(seq2)] \
                      if tm["type"]=="subpop_ell" else g_cols[gi%len(g_cols)]
                cscale[ti] = [[0,col],[1,col]]

        relayout = {
            "paper_bgcolor":p["bg"],"plot_bgcolor":p["bg"],
            "font.color":p["text"],"title.font.color":p["text"],
            "scene.bgcolor":p["panel"],
            **{f"scene.{ax}.backgroundcolor":p["panel"]
               for ax in ["xaxis","yaxis","zaxis"]},
            **{f"scene.{ax}.gridcolor":p["grid"]
               for ax in ["xaxis","yaxis","zaxis"]},
            **{f"scene.{ax}.zerolinecolor":p["grid"]
               for ax in ["xaxis","yaxis","zaxis"]},
            **{f"scene.{ax}.tickfont.color":p["subtext"]
               for ax in ["xaxis","yaxis","zaxis"]},
            **{f"scene.{ax}.title.font.color":p["text"]
               for ax in ["xaxis","yaxis","zaxis"]},
            "legend.bgcolor":p["bg"],"legend.bordercolor":p["grid"],
            "legend.font.color":p["text"],"legend.title.font.color":p["text"],
        }
        palette_buttons.append(dict(
            label=pal_name, method="update",
            args=[{"marker.color":safe_c,"colorscale":cscale},relayout],
        ))

    active_idx = list(PALETTES.keys()).index(
        f"{pal['mode']} / {pal['background']}"
    )

    # ── Layout ───────────────────────────────────────────────────────
    axis_style = dict(
        backgroundcolor=pal["panel"],gridcolor=pal["grid"],
        showbackground=True,zerolinecolor=pal["grid"],
        tickfont=dict(color=pal["subtext"],size=9),
    )
    fig.update_layout(
        paper_bgcolor=pal["bg"],plot_bgcolor=pal["bg"],
        font=dict(color=pal["text"]),
        title=dict(
            text=(f"<b>CRISP Step 6 — PCA Pass 1 — 3D Scatter</b><br>"
                  f"<sup>N={n_tot} — {n_grp} group(s) — {n_out} outlier(s) "
                  f"— Confidence: PROBABLE</sup>"),
            x=0.0,xanchor="left",font=dict(color=pal["text"],size=13),
        ),
        scene=dict(
            bgcolor=pal["panel"],
            xaxis=dict(title=dict(text=f"PC{args.pc_x}",
                       font=dict(color=pal["text"],size=10)),**axis_style),
            yaxis=dict(title=dict(text=f"PC{args.pc_y}",
                       font=dict(color=pal["text"],size=10)),**axis_style),
            zaxis=dict(title=dict(text=f"PC{args.pc_z}",
                       font=dict(color=pal["text"],size=10)),**axis_style),
            camera=dict(eye=dict(x=1.5,y=1.5,z=0.8)),
        ),
        legend=dict(
            x=1.01,y=1.0,xanchor="left",yanchor="top",
            bgcolor=pal["bg"],bordercolor=pal["grid"],borderwidth=1,
            font=dict(color=pal["text"],size=9),
            title=dict(text="Group",font=dict(color=pal["text"],size=10)),
        ),
        margin=dict(l=0,r=220,t=80,b=20),
        width=1100,height=750,
        updatemenus=[
            dict(type="buttons",direction="right",
                 x=0.0,xanchor="left",y=-0.04,yanchor="top",
                 showactive=True,active=active_idx,
                 bgcolor=pal["panel"],bordercolor=pal["grid"],
                 font=dict(color=pal["text"],size=10),
                 buttons=palette_buttons,pad=dict(r=4,t=0,b=0)),
        ],
        annotations=[
            dict(text="<b>Palette:</b>",xref="paper",yref="paper",
                 x=0,y=-0.04,showarrow=False,
                 font=dict(size=10,color=pal["subtext"]),
                 align="left",xanchor="right"),
        ],
    )

    return fig, trace_meta, all_groups, subpop_map, centroids


##########################################################################
# JS INJECTION
##########################################################################

def inject_js(html, trace_meta, all_groups, subpop_map,
              centroids, df, pal, args):
    import re as _re

    # WHY regex injection:
    # Plotly's write_html silently drops visible=False for Scatter3d traces.
    # We post-process the HTML JSON directly to inject "visible":false
    # for any trace whose name contains special unicode chars (›, ◆, ▲)
    # that identify sub-pop and centroid traces.
    def _inject_visible_false(h):
        pattern = r'("name":"[^"]*(?:\\u203a|\\u25c6|\\u25b2|\\u25b3|\\u2ca)[^"]*")(,")' 
        def repl(m):
            return m.group(1) + ',"visible":false' + m.group(2)
        return _re.sub(pattern, repl, h)
    html = _inject_visible_false(html)

    pc_x = f"PC{args.pc_x}"
    pc_y = f"PC{args.pc_y}"
    pc_z = f"PC{args.pc_z}"

    # Per-group bounding boxes
    grp_ranges = {}
    for grp in all_groups:
        gd = df[(df["SUPERPOP"]==grp)&(df["OUTLIER"]!="YES")]
        if len(gd)==0: continue
        pad = 0.15
        grp_ranges[grp] = {
            "x":[float(gd[pc_x].min()*(1+pad)),float(gd[pc_x].max()*(1+pad))],
            "y":[float(gd[pc_y].min()*(1+pad)),float(gd[pc_y].max()*(1+pad))],
            "z":[float(gd[pc_z].min()*(1+pad)),float(gd[pc_z].max()*(1+pad))],
        }

    # Per-subpop bounding boxes
    sub_ranges = {}
    for grp,subs in subpop_map.items():
        sub_ranges[grp] = {}
        for sub in subs:
            sd = df[(df["SUPERPOP"]==grp)&(df["SUBPOP"]==sub)&
                    (df["OUTLIER"]!="YES")]
            if len(sd)==0: continue
            pad = 0.20
            sub_ranges[grp][sub] = {
                "x":[float(sd[pc_x].min()*(1+pad)),float(sd[pc_x].max()*(1+pad))],
                "y":[float(sd[pc_y].min()*(1+pad)),float(sd[pc_y].max()*(1+pad))],
                "z":[float(sd[pc_z].min()*(1+pad)),float(sd[pc_z].max()*(1+pad))],
            }

    js_data = json.dumps({
        "trace_meta" : trace_meta,
        "all_groups" : all_groups,
        "subpop_map" : subpop_map,
        "grp_ranges" : grp_ranges,
        "sub_ranges" : sub_ranges,
        "pal"        : {k:pal[k] for k in
                        ["bg","panel","grid","text","subtext","fail"]},
    })

    js_block = f"""
<style>
#crisp-ctrl {{
    font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;
    padding:10px 16px 16px 16px;
    background:{pal['bg']};
    border-top:1px solid {pal['grid']};
    user-select:none;
}}
.ctrl-section {{ margin-bottom:10px; }}
.ctrl-row {{ display:flex; flex-wrap:wrap; align-items:center; gap:5px; margin-bottom:5px; }}
.ctrl-label {{ font-size:11px; font-weight:600; color:{pal['subtext']}; min-width:130px; }}
.crisp-btn {{
    font-size:10px; padding:4px 10px; border-radius:5px;
    border:1px solid {pal['grid']}; cursor:pointer;
    background:{pal['panel']}; color:{pal['text']};
    transition:all 0.15s ease;
}}
.crisp-btn:hover {{ background:{pal['grid']}; }}
.crisp-btn.active {{ background:{pal['text']}; color:{pal['bg']}; border-color:{pal['text']}; font-weight:600; }}
.crisp-btn.action {{ background:{pal['panel']}; }}
.crisp-btn.danger {{ background:{pal['fail']}22; border:1px solid {pal['fail']}88; color:{pal['text']}; }}
.crisp-btn.focus-btn {{ border:1px dashed {pal['grid']}; }}
.crisp-btn.focus-btn.active {{ background:{pal['text']}22; border:1px dashed {pal['text']}; font-weight:600; }}
#sub-panel {{ display:none; margin-top:6px; padding:8px 10px;
    border:1px solid {pal['grid']}; border-radius:6px; background:{pal['panel']}22; }}
.sub-mode-note {{ font-size:10px; color:{pal['subtext']}; font-style:italic; margin-bottom:4px; }}
</style>

<div id="crisp-ctrl">
  <div class="ctrl-section">
    <div class="ctrl-row">
      <span class="ctrl-label">Groups (multi-select):</span>
      <button class="crisp-btn action" title="Show all population groups" onclick="crispShowAll()">Show All</button>
      <button class="crisp-btn action" title="Deselect all groups" onclick="crispClear()">Clear All</button>
      <button class="crisp-btn danger" title="Zoom axes to bounding box of selected groups" onclick="crispRescale()">⟳ Rescale to selection</button>
      <button class="crisp-btn action" title="Hide all groups not in current selection" onclick="crispIsolateSelection()">◎ Isolate selection</button>
      <button class="crisp-btn action" title="Reset axes, camera and restore all groups" onclick="crispResetView()">↺ Reset view</button>
    </div>
    <div class="ctrl-row" id="grp-btn-row"></div>
  </div>
  <div class="ctrl-section">
    <div class="ctrl-row">
      <span class="ctrl-label">Focus (isolate + zoom):</span>
      <button class="crisp-btn action" title="Reset — restore all groups" onclick="crispResetFocus()">↺ Reset</button>
      <span id="focus-btn-row" style="display:contents"></span>
    </div>
  </div>
  <div id="sub-panel">
    <div class="sub-mode-note" id="sub-mode-note"></div>
    <div class="ctrl-row">
      <span class="ctrl-label">Sub-populations:</span>
      <button class="crisp-btn action" title="Show all sub-populations" onclick="crispSubShowAll()">Show All</button>
      <button class="crisp-btn action" title="Hide all sub-populations" onclick="crispSubClear()">Clear</button>
      <button class="crisp-btn danger" title="Zoom axes to selected sub-populations" onclick="crispSubRescale()">⟳ Rescale to sub-selection</button>
    </div>
    <div class="ctrl-row" id="sub-btn-row"></div>
  </div>
</div>

<script>
(function() {{
    const allGroups  = {json.dumps(all_groups)};
    const subpopMap  = {json.dumps(subpop_map)};
    const grpRanges  = {json.dumps(grp_ranges)};
    const subRanges  = {json.dumps(sub_ranges)};
    const plotDiv    = document.querySelector('.plotly-graph-div');

    let selectedGroups = new Set(allGroups);
    let selectedSubs   = new Set();
    let focusedGroup   = null;

    // ── Name-based trace classifier ──────────────────────────────────
    // WHY name-based: legendgroup values are misaligned in Plotly's
    // serialised JSON output. Trace names are correctly aligned.
    // All classification uses trace name patterns only.
    function classifyTrace(name) {{
        const n = name || '';
        if (n.startsWith('Outlier'))                          return {{type:'outlier',  grp:'',    sub:''}};
        if (n.includes(' ▲'))                                 return {{type:'centroid', grp:n.replace(' ▲',''), sub:''}};
        if (n.includes(' › ') && n.includes('(n='))          return {{type:'sub_col',  grp:n.split(' › ')[0], sub:n.split(' › ')[1].split(' (')[0]}};
        if (n.includes(' › ') && n.includes('◆'))            return {{type:'sub_shp',  grp:n.split(' › ')[0], sub:n.split(' › ')[1].replace(' ◆','').trim()}};
        if (n.endsWith('_ell') && n.includes('_') && n.split('_').length>2) {{
            const parts = n.replace('_ell','').split('_');
            return {{type:'sub_ell', grp:parts[0], sub:parts.slice(1).join('_')}};
        }}
        if (n.endsWith('_ell'))                               return {{type:'grp_ell', grp:n.replace('_ell',''), sub:''}};
        if (n.includes('(n='))                                return {{type:'scatter',  grp:n.split(' (n=')[0], sub:''}};
        return {{type:'unknown', grp:'', sub:''}};
    }}

    // ── Core visibility setter ────────────────────────────────────────
    function setVis(visFn, showLegFn) {{
        const n = plotDiv.data.length;
        const vis     = [];
        const showLeg = [];
        for (let i=0;i<n;i++) {{
            const t  = plotDiv.data[i];
            const cl = classifyTrace(t.name);
            vis.push(visFn(cl, t.name));
            showLeg.push(showLegFn ? showLegFn(cl, t.name) : undefined);
        }}
        // WHY no extra array wrapping:
        // Plotly.restyle expects visible as a flat boolean array.
        // Wrapping in an extra [] creates a nested array Plotly ignores.
        // This was the root cause of Focus not working.
        const update = {{'visible': vis}};
        if (showLegFn) update['showlegend'] = showLeg;
        Plotly.restyle(plotDiv, update);
    }}

    // ── Render group buttons ─────────────────────────────────────────
    function renderGroupButtons() {{
        const row = document.getElementById('grp-btn-row');
        row.innerHTML = '';
        allGroups.forEach(grp => {{
            const btn = document.createElement('button');
            const on  = selectedGroups.has(grp);
            btn.className = 'crisp-btn' + (on ? ' active' : '');
            btn.textContent = grp;
            btn.title = on ? `Deselect ${{grp}}` : `Add ${{grp}} to selection`;
            btn.onclick = () => {{ toggleGroup(grp); }};
            row.appendChild(btn);
        }});
    }}

    function renderFocusButtons() {{
        const row = document.getElementById('focus-btn-row');
        row.innerHTML = '';
        allGroups.forEach(grp => {{
            const btn = document.createElement('button');
            const on  = focusedGroup === grp;
            btn.className = 'crisp-btn focus-btn' + (on ? ' active' : '');
            btn.textContent = `⊙ ${{grp.slice(0,11)}}`;
            btn.title = `Isolate ${{grp}} and zoom — click again to reset`;
            btn.onclick = () => {{
                if (focusedGroup === grp) {{ crispResetFocus(); }}
                else {{ focusGroup(grp); }}
            }};
            row.appendChild(btn);
        }});
    }}

    function toggleGroup(grp) {{
        if (selectedGroups.has(grp)) {{ selectedGroups.delete(grp); }}
        else {{ selectedGroups.add(grp); }}
        applyGroupVisibility();
        updateSubPanel();
        renderGroupButtons();
    }}

    function applyGroupVisibility() {{
        setVis(
            cl => {{
                if (cl.type==='outlier')  return false; // hide outliers in normal view
                if (cl.type==='sub_col' || cl.type==='sub_shp' || cl.type==='sub_ell') return false;
                if (cl.type==='centroid') return selectedGroups.has(cl.grp);
                if (cl.type==='grp_ell') return false; // hide ellipsoids in normal view
                if (cl.type==='scatter') return selectedGroups.has(cl.grp);
                return false;
            }},
            cl => cl.type==='scatter'  // only scatter traces get legend entries
        );
    }}

    window.crispShowAll = function() {{
        selectedGroups = new Set(allGroups);
        applyGroupVisibility();
        updateSubPanel();
        renderGroupButtons();
    }};

    window.crispClear = function() {{
        selectedGroups = new Set();
        applyGroupVisibility();
        updateSubPanel();
        renderGroupButtons();
    }};

    window.crispRescale = function() {{
        if (selectedGroups.size===0) return;
        let xMin=Infinity,xMax=-Infinity,yMin=Infinity,yMax=-Infinity,zMin=Infinity,zMax=-Infinity;
        selectedGroups.forEach(g => {{
            const r=grpRanges[g]; if(!r) return;
            xMin=Math.min(xMin,r.x[0]);xMax=Math.max(xMax,r.x[1]);
            yMin=Math.min(yMin,r.y[0]);yMax=Math.max(yMax,r.y[1]);
            zMin=Math.min(zMin,r.z[0]);zMax=Math.max(zMax,r.z[1]);
        }});
        Plotly.relayout(plotDiv,{{'scene.xaxis.range':[xMin,xMax],'scene.yaxis.range':[yMin,yMax],'scene.zaxis.range':[zMin,zMax]}});
    }};

    window.crispIsolateSelection = function() {{
        if (selectedGroups.size===0) return;
        setVis(
            cl => {{
                if (cl.type==='sub_col'||cl.type==='sub_shp'||cl.type==='sub_ell') return false;
                if (cl.type==='grp_ell') return false;
                if (cl.type==='outlier') return false;
                return selectedGroups.has(cl.grp);
            }},
            cl => cl.type==='scatter' && selectedGroups.has(cl.grp)
        );
    }};

    // ── Focus: isolate single group + zoom ──────────────────────────
    function focusGroup(grp) {{
        focusedGroup = grp;
        renderFocusButtons();

        // Show ONLY this group's scatter trace
        // Sub-pop traces handled by updateSubPanel below
        setVis(
            cl => {{
                if (cl.type==='scatter' && cl.grp===grp) return true;
                return false;
            }},
            cl => cl.type==='scatter' && cl.grp===grp
        );

        // Zoom camera + axes
        const r = grpRanges[grp];
        const cent = {json.dumps({g: list(v) for g,v in centroids.items()})};
        const relayout = {{}};
        if (cent[grp]) {{
            const [cx,cy,cz] = cent[grp];
            relayout['scene.camera'] = {{center:{{x:cx*50,y:cy*50,z:cz*50}},eye:{{x:0.8,y:0.8,z:0.5}}}};
        }}
        if (r) {{
            relayout['scene.xaxis.range']=r.x;
            relayout['scene.yaxis.range']=r.y;
            relayout['scene.zaxis.range']=r.z;
        }}
        if (Object.keys(relayout).length>0) Plotly.relayout(plotDiv, relayout);

        // Select this group and open sub-panel
        selectedGroups = new Set([grp]);
        selectedSubs   = new Set();
        renderGroupButtons();
        updateSubPanel();
    }}

    window.crispResetFocus = function() {{
        focusedGroup   = null;
        selectedGroups = new Set(allGroups);
        selectedSubs   = new Set();
        applyGroupVisibility();
        Plotly.relayout(plotDiv,{{
            'scene.camera':{{eye:{{x:1.5,y:1.5,z:0.8}},center:{{x:0,y:0,z:0}}}},
            'scene.xaxis.range':undefined,'scene.yaxis.range':undefined,'scene.zaxis.range':undefined,
        }});
        document.getElementById('sub-panel').style.display='none';
        renderGroupButtons();
        renderFocusButtons();
    }};

    window.crispResetView = window.crispResetFocus;

    // ── Sub-population panel ─────────────────────────────────────────
    function updateSubPanel() {{
        const panel = document.getElementById('sub-panel');
        const grpsWithSubs = Array.from(selectedGroups).filter(g=>subpopMap[g]&&subpopMap[g].length>0);
        if (grpsWithSubs.length===0) {{
            panel.style.display='none';
            return;
        }}
        panel.style.display='block';
        const isSingle = selectedGroups.size===1;
        document.getElementById('sub-mode-note').textContent = isSingle
            ? 'Single group — sub-populations shown in distinct colours'
            : 'Multiple groups — sub-populations shown as distinct shapes within group colour';

        const allSubs = grpsWithSubs.flatMap(g=>subpopMap[g].map(s=>({{grp:g,sub:s}})));
        if (selectedSubs.size===0) selectedSubs = new Set(allSubs.map(x=>x.grp+'|'+x.sub));
        renderSubButtons(allSubs, isSingle);
        applySubVisibility(allSubs, isSingle);
    }}

    function renderSubButtons(allSubs, isSingle) {{
        const row = document.getElementById('sub-btn-row');
        row.innerHTML='';
        allSubs.forEach(({{grp,sub}}) => {{
            const key = grp+'|'+sub;
            const on  = selectedSubs.has(key);
            const btn = document.createElement('button');
            btn.className = 'crisp-btn'+(on?' active':'');
            btn.textContent = sub;
            btn.title = `Toggle ${{grp}} › ${{sub}} — ${{isSingle?'distinct colour':'distinct shape'}}`;
            btn.onclick = () => {{
                if (selectedSubs.has(key)) selectedSubs.delete(key);
                else selectedSubs.add(key);
                renderSubButtons(allSubs, isSingle);
                applySubVisibility(allSubs, isSingle);
            }};
            row.appendChild(btn);
        }});
    }}

    function applySubVisibility(allSubs, isSingle) {{
        const activeKeys = new Set(selectedSubs);
        setVis(
            cl => {{
                // Group scatter: hide when sub-pops are active for it
                if (cl.type==='scatter') {{
                    if (!selectedGroups.has(cl.grp)) return false;
                    const hasSubs = subpopMap[cl.grp]&&subpopMap[cl.grp].length>0;
                    return !(hasSubs && activeKeys.size>0);
                }}
                // Centroid: show for focused group
                if (cl.type==='centroid') return selectedGroups.has(cl.grp);
                // Group ellipsoid: hide when sub-pops active
                if (cl.type==='grp_ell') return false;
                // Sub colour: single group mode
                if (cl.type==='sub_col') {{
                    return isSingle && selectedGroups.has(cl.grp) && activeKeys.has(cl.grp+'|'+cl.sub);
                }}
                // Sub shape: multi group mode
                if (cl.type==='sub_shp') {{
                    return !isSingle && selectedGroups.has(cl.grp) && activeKeys.has(cl.grp+'|'+cl.sub);
                }}
                // Sub ellipsoid: single group mode only
                if (cl.type==='sub_ell') {{
                    return isSingle && selectedGroups.has(cl.grp) && activeKeys.has(cl.grp+'|'+cl.sub);
                }}
                if (cl.type==='outlier') return false;
                return false;
            }},
            cl => {{
                if (cl.type==='sub_col') return isSingle && selectedGroups.has(cl.grp) && activeKeys.has(cl.grp+'|'+cl.sub);
                if (cl.type==='sub_shp') return !isSingle && selectedGroups.has(cl.grp) && activeKeys.has(cl.grp+'|'+cl.sub);
                if (cl.type==='scatter') return false; // hide group scatter in legend when sub-pops shown
                return false;
            }}
        );
    }}

    window.crispSubShowAll = function() {{
        const grpsWithSubs = Array.from(selectedGroups).filter(g=>subpopMap[g]);
        grpsWithSubs.forEach(g=>subpopMap[g].forEach(s=>selectedSubs.add(g+'|'+s)));
        const isSingle = selectedGroups.size===1;
        const allSubs  = grpsWithSubs.flatMap(g=>subpopMap[g].map(s=>({{grp:g,sub:s}})));
        renderSubButtons(allSubs, isSingle);
        applySubVisibility(allSubs, isSingle);
    }};

    window.crispSubClear = function() {{
        selectedSubs = new Set();
        const grpsWithSubs = Array.from(selectedGroups).filter(g=>subpopMap[g]);
        const isSingle = selectedGroups.size===1;
        const allSubs  = grpsWithSubs.flatMap(g=>subpopMap[g].map(s=>({{grp:g,sub:s}})));
        renderSubButtons(allSubs, isSingle);
        applySubVisibility(allSubs, isSingle);
    }};

    window.crispSubRescale = function() {{
        if (selectedSubs.size===0) return;
        let xMin=Infinity,xMax=-Infinity,yMin=Infinity,yMax=-Infinity,zMin=Infinity,zMax=-Infinity;
        selectedSubs.forEach(key => {{
            const [grp,sub]=key.split('|');
            const r=(subRanges[grp]||{{}})[sub]; if(!r) return;
            xMin=Math.min(xMin,r.x[0]);xMax=Math.max(xMax,r.x[1]);
            yMin=Math.min(yMin,r.y[0]);yMax=Math.max(yMax,r.y[1]);
            zMin=Math.min(zMin,r.z[0]);zMax=Math.max(zMax,r.z[1]);
        }});
        if (xMin===Infinity) return;
        Plotly.relayout(plotDiv,{{'scene.xaxis.range':[xMin,xMax],'scene.yaxis.range':[yMin,yMax],'scene.zaxis.range':[zMin,zMax]}});
    }};

    // ── Init ─────────────────────────────────────────────────────────
    renderGroupButtons();
    renderFocusButtons();
    // Start with clean state — no ellipsoids, no sub-pops
    applyGroupVisibility();
}})();
</script>
"""

    return html.replace("</body>", js_block + "\n</body>")


##########################################################################
# STATIC BACKUP
##########################################################################

def plot_static(df, pal, args, all_groups, grp_cols, out_path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa
    except ImportError:
        return
    pc_x,pc_y,pc_z = f"PC{args.pc_x}",f"PC{args.pc_y}",f"PC{args.pc_z}"
    views = [(30,45,"View 1"),(20,135,"View 2"),(60,90,"View 3")]
    fig   = plt.figure(figsize=(18,6),facecolor=pal["bg"])
    for idx,(elev,azim,label) in enumerate(views,1):
        ax = fig.add_subplot(1,3,idx,projection="3d",facecolor=pal["panel"])
        for grp in all_groups:
            gd = df[(df["SUPERPOP"]==grp)&(df["OUTLIER"]!="YES")]
            if len(gd)==0: continue
            ax.scatter(gd[pc_x],gd[pc_y],gd[pc_z],
                       c=grp_cols[grp],s=3,alpha=0.55,
                       label=grp if idx==1 else None)
        out = df[df["OUTLIER"]=="YES"]
        if len(out)>0:
            ax.scatter(out[pc_x],out[pc_y],out[pc_z],
                       c=pal["fail"],s=18,alpha=0.9,marker="x",
                       linewidths=1.2,
                       label="Outlier" if idx==1 else None)
        ax.set_xlabel(pc_x,color=pal["text"],fontsize=8)
        ax.set_ylabel(pc_y,color=pal["text"],fontsize=8)
        ax.set_zlabel(pc_z,color=pal["text"],fontsize=8)
        ax.set_title(label,color=pal["text"],fontsize=9)
        ax.view_init(elev=elev,azim=azim)
        ax.tick_params(colors=pal["subtext"],labelsize=6)
        ax.xaxis.pane.fill=False
        ax.yaxis.pane.fill=False
        ax.zaxis.pane.fill=False
    h,l = fig.axes[0].get_legend_handles_labels()
    fig.legend(h,l,loc="center right",bbox_to_anchor=(1.0,0.5),
               frameon=True,facecolor=pal["bg"],edgecolor=pal["grid"],
               fontsize=7,title="Group",title_fontsize=9)
    fig.suptitle(f"CRISP Step 6 — PCA 3D Static — "
                 f"PC{args.pc_x}×PC{args.pc_y}×PC{args.pc_z}",
                 color=pal["text"],fontsize=10,y=1.02)
    fig.savefig(out_path,bbox_inches="tight",dpi=200,facecolor=pal["bg"])
    plt.close(fig)
    print(f"[CRISP] Static PDF: {out_path}")


##########################################################################
# MAIN
##########################################################################

def main():
    args    = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[CRISP plot_pca_pass1_3d.py] v1.4.0")
    print(f"[CRISP plot_pca_pass1_3d.py] Project : {args.project}")
    print(f"[CRISP plot_pca_pass1_3d.py] Axes    : "
          f"PC{args.pc_x} × PC{args.pc_y} × PC{args.pc_z}")

    pal = active_palette()
    print(f"[CRISP plot_pca_pass1_3d.py] Palette : "
          f"{pal['mode']} / {pal['background']}")

    df = load_data(args)
    all_groups = sorted(df["SUPERPOP"].unique().tolist())
    seq = pal["group_cols"]
    grp_cols = {g:seq[i%len(seq)] for i,g in enumerate(all_groups)}

    fig, trace_meta, all_groups, subpop_map, centroids = build_figure(df, args)

    # WHY explicit update_traces:
    # Plotly's write_html silently drops visible=False for some trace types
    # (Scatter3d with visible=False is particularly affected).
    # We force-set visible on every trace by name pattern before writing HTML.
    for i, trace in enumerate(fig.data):
        name = trace.name or ""
        # Hide: sub-pop variants, ellipsoids, centroids, outliers
        if ("›" in name or "_ell" in name or
            "◆" in name or " ▲" in name or
            name.startswith("Outlier")):
            fig.data[i].visible = False
        else:
            # Group scatter traces start visible
            fig.data[i].visible = True

    out_html = out_dir / f"{args.project}_pca_3d.html"
    fig.write_html(
        str(out_html), include_plotlyjs="cdn", full_html=True,
        config=dict(displayModeBar=True,
                    modeBarButtonsToRemove=["toImage"],
                    displaylogo=False, scrollZoom=True),
    )

    # Inject JS
    with open(out_html) as f:
        html = f.read()
    html = inject_js(html, trace_meta, all_groups,
                     subpop_map, centroids, df, pal, args)
    with open(out_html,"w") as f:
        f.write(html)

    print(f"[CRISP plot_pca_pass1_3d.py] HTML: {out_html}")

    if args.static.upper()=="YES":
        out_pdf = out_dir / f"{args.project}_pca_3d_static.pdf"
        plot_static(df, pal, args, all_groups, grp_cols, str(out_pdf))

    print(f"[CRISP plot_pca_pass1_3d.py] Complete.")


if __name__ == "__main__":
    main()
