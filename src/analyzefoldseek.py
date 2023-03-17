#!/usr/bin/env python3

import pandas as pd
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import sys

def getstats(df, metrics = ['evalue', "alnlen", "pident", "alntmscore"]):
    return (
        df[metrics]
        .agg(["count", np.mean, np.std, np.min, np.median, np.max])
        .transpose()
    )

def main():
    print("analyzefoldseek.py")
    # Load results database and preprocess
    aln = pd.read_csv(
        sys.argv[1],
        sep='\t',
        names = "query,target,pident,alnlen,evalue,bits,qlen,tlen,qstart,qend,tstart,tend,qaln,taln,qseq,tseq,alntmscore".split(','),
    )
    aln["query"] = aln["query"].str.strip(".pdb").str.upper()
    aln["target"] = aln["target"].str.strip(".pdb").str.upper()

    print()
    print("Table. Hits statistic for relevant metrics")
    print(getstats(aln).to_markdown(floatfmt=".3g"))
    print()

    # Plot results
    aln["p-evalue"] = -np.log10(aln["evalue"])
    fig = make_subplots(
    rows=1, cols=3)

    fig.add_trace(go.Scatter(
        x = aln["alnlen"],
        y = aln["pident"],
        mode='markers',
        customdata = aln["target"],
        hovertemplate="<br>Target:%{customdata}",
        marker_color=aln["p-evalue"], 
    ),
    row=1, col=1
    )

    fig.add_trace(go.Scatter(
        x = aln["alnlen"],
        y = aln["alntmscore"],
        mode='markers',
        customdata = aln["target"],
        hovertemplate="<br>Target:%{customdata}",
        marker_color=aln["p-evalue"],
    ),
    row=1, col=2
    )
    fig.add_trace(go.Scatter(
        x = aln["pident"],
        y = aln["alntmscore"],
        mode='markers',
        customdata = aln["target"],
        hovertemplate="<br>Target:%{customdata}",
        marker=dict(colorbar=dict(title = "p-evalue"), color=aln["p-evalue"]),
    ),
    row=1, col=3
    )

    # Update xaxis properties
    fig.update_xaxes(title_text="Alignment Lenght", row=1, col=1)
    fig.update_xaxes(title_text="Alignment Lenght", row=1, col=2)
    fig.update_xaxes(title_text="Percent Identity", row=1, col=3)

    # Update yaxis properties
    fig.update_yaxes(title_text="Percent Identity", row=1, col=1)
    fig.update_yaxes(title_text="Alignment TMscore", row=1, col=2)
    fig.update_yaxes(title_text="Alignment TMscore", row=1, col=3)

    fig.update_layout(height=600, width=1700,
                  title_text="Foldseek analysis", showlegend=False)
    fig.write_html("Regressions.html")

main()
