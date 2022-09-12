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