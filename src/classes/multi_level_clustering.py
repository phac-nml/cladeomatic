import pandas as pd
import scipy

class multi_level_clustering:
    cluster_memberships = {}
    thresholds = []
    labels = []
    linkage = None

    def __init__(self,dist_mat_file,thresholds,method):
        df = self.read_matrix(dist_mat_file).astype(float)
        self.labels = df.columns.values.tolist()
        matrix = scipy.spatial.distance.squareform(df.values)
        del df
        self.thresholds = thresholds
        self.linkage = scipy.cluster.hierarchy.linkage(matrix, method=method, metric='precomputed')
        self.init_membership()
        self.assign_clusters()

    def init_membership(self):
        for label in self.labels:
            self.cluster_memberships[label] = []

    def read_matrix(self,f):
        return pd.read_csv(f,header=0,sep="\t",index_col=0)

    def assign_clusters(self):
        for idx,dist in enumerate(self.thresholds):
            cid = 0
            clusters = scipy.cluster.hierarchy.fcluster(self.linkage, dist, criterion='distance')
            for label in self.labels:
                self.cluster_memberships[label].append(f'{clusters[cid]}')
                cid+=1


    def get_memberships(self):
        return self.cluster_memberships


