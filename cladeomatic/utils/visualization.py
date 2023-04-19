import sys

import plotly.express as px
import pandas as pd
from deprecated import deprecated


@deprecated()
def plot_bar(x,y):
    """
    A method to plot a bar chart with the x and y objects passed
    :param x: obj - list or a dictionary for the x axis objects
    :param y: obj - list or a dictionary for the y axis objects
    :return: Figure - return the figure produced
    """
    if len(y) == 0:
        return None
    #create the pandas data frome from a zip of the x and y objects
    #cast to a list
    df = pd.DataFrame(list(zip(x,y)),columns=['x','y'])
    fig = px.bar(df, x="x",
                 y='y',

                 )
    fig.update_layout(yaxis_range=[0, max(y)])
    return fig

def create_dist_histo(data,outfile):
    """
    A method to create a histogram from the data dictionary passed.
    This method also writes the histogram to a file.
    :param data: dictionary - the data to be plotted
    :param outfile: String - the output figure file path
    """
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
