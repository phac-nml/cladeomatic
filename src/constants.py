from src.version import __version__
EXTENSIONS = {'text': ['.txt','.tsv','.mat','.text'],
    'parquet': ['.parq','.parquet','.pq']}
PD_HEADER = [
    'query_id',
    'ref_id',
    'dist'
]

MIN_FILE_SIZE = 32


CLUSTER_METHODS = ['average','complete','single']


MC_RUN_DATA = {
    'genomic address service: de novo clustering': f'version: {__version__}',
    'analysis_start_time':'',
    'analysis_end_time':'',
    'parameters':{},
    'threshold_map':{},
    'result_file':''
}
