import sys

import plotly.express as px
import pandas as pd


def plot_bar(x,y):
    if len(y) == 0:
        return None
    df = pd.DataFrame(list(zip(x,y)),columns=['x','y'])
    fig = px.bar(df, x="x",
                 y='y',

                 )
    fig.update_layout(yaxis_range=[0, max(y)])
    return fig

def create_dist_histo(data,outfile):
    fig = px.bar(x=data.keys(),y=data.values())
    fig.update_yaxes(title="Frequency")
    fig.update_xaxes(title="Sample Hamming SNP distance")
    fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',
        template="simple_white",
    )
    fh = open(outfile, 'w')
    fh.write(fig.to_html())
    fh.close()
