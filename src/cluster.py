import os.path
import sys
import json
from datetime import datetime
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from src.version import __version__
from src.constants import CLUSTER_METHODS, MC_RUN_DATA
from src.classes.multi_level_clustering import multi_level_clustering
from src.utils import is_file_ok

def parse_args():
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass

    parser = ArgumentParser(
        description="Genomic Address Service: De novo hierarchical sequence clustering",
        formatter_class=CustomFormatter)
    parser.add_argument('-i','--matrix', type=str, required=True,help='TSV formated distance matrix or parquet')
    parser.add_argument('-o','--outdir', type=str, required=True, help='Output directory to put cluster results')
    parser.add_argument('-m','--method', type=str, required=False, help='cluster method [single, complete, average]',default='average')
    parser.add_argument('-t','--thresholds', type=str, required=True, help='thresholds delimited by ,')
    parser.add_argument('-d', '--delimeter', type=str, required=False, help='delimeter desired for nomenclature code',default=".")
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    parser.add_argument('-f', '--force', required=False, help='Overwrite existing directory',
                        action='store_true')

    return parser.parse_args()


def write_clusters(clusters,num_thresholds,file,delimeter="."):
    header = ['id','address']
    for i in range(num_thresholds):
        header.append(f'level_{i+1}')
    header = "\t".join([str(x) for x in header])
    with open(file,'w') as fh:
        fh.write("{}\n".format(header))
        for id in clusters:
            address = f'{delimeter}'.join([str(x) for x in clusters[id]])
            fh.write("{}\n".format("\t".join(str(x) for x in ([id, address ] + clusters[id]))))

def format_threshold_map(thresholds):
    data = {}
    for i,value in enumerate(thresholds):
        data[f'level_{i+1}'] = value
    return data

def write_threshold_map(data,file):
    with open(file,'w') as fh:
        fh.write(json.dumps(data, indent=4))

    fh.close()


def run():
    cmd_args = parse_args()
    matrix = cmd_args.matrix
    outdir = cmd_args.outdir
    method = cmd_args.method
    thresholds = [float(x) for x in cmd_args.thresholds.split(',')]
    delimeter= cmd_args.delimeter
    force = cmd_args.force

    run_data = MC_RUN_DATA
    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = vars(cmd_args)
    t_map = format_threshold_map(thresholds)
    run_data['threshold_map'] = t_map
    if not is_file_ok(matrix):
        print(f'Error {matrix} does not exist or is empty')
        sys.exit()

    if not method in CLUSTER_METHODS:
        print(f'Error {method} is not one of the accepeted methods {CLUSTER_METHODS}')
        sys.exit()

    if os.path.isdir(outdir) and not force:
        print(f'Error {outdir} exists, if you would like to overwrite, then specify --force')
        sys.exit()

    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)

    mc = multi_level_clustering(matrix,thresholds,method)

    memberships = mc.get_memberships()

    if len(memberships) == 0:
        print(f'Error something when wrong during clustering')
        sys.exit()

    run_data['result_file'] = os.path.join(outdir,"clusters.text")

    write_clusters(memberships, len(thresholds), run_data['result_file'], delimeter)

    write_threshold_map(t_map, os.path.join(outdir,"thresholds.json"))

    run_data['analysis_end_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    with open(os.path.join(outdir,"run.json"),'w') as fh:
        fh.write(json.dumps(run_data, indent=4))


# call main function
if __name__ == '__main__':
    run()
