import plotly.express as px

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

