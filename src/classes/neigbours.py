import os
from numba import jit
from numba.typed import List
import pyarrow.parquet as pq

class neighbours:
    file_path = None
    threshhold = 0
    batch_size = None
    num_lines = 0
    delim = "\t"
    is_ok = True
    error_msg = []
    components = []
    sample_ids = []
    near_neighbors = {}

    def __init__(self,file_path,out_path,is_square=True,thresh=0,delim='\t'):
        self.file_path = file_path
        self.out_path = out_path
        self.threshhold = thresh
        self.delim = delim

        if not os.path.isfile(file_path):
            self.error_msg.append("Error matrix file: {} does not exist".format(self.file_path))
            self.is_ok = False
            return None

        self.num_lines = self.get_file_length() - 1

        if self.num_lines < 2:
            self.error_msg.append("Error matrix file: {} does not contain at least two samples".format(self.file_path))
            self.is_ok = False
            return None

        if is_square:
            self.parse_square_distance_matrix()
        else:
            self.parse_ltriangle_distance_matrix()

        c = len(self.components)

        self.link()
        self.link()
        while c != len(self.components):
            c = len(self.components)
            self.link()

        self.write_links()

        return None


    def get_result(self):
        return {
            'threshold':self.threshhold,
            'components':self.components,
            'sample_ids':self.sample_ids,
            'near_neighbors':self.near_neighbors
        }

    def get_file_length(self):
        return int(os.popen(f'wc -l {self.file_path}').read().split()[0])

    def init_nearest_neighbor(self,samples):
        for sample_id in samples:
            self.near_neighbors[sample_id] = {
                'dist':None,
                'sample_ids':[]
            }

    def parse_square_distance_matrix(self):
        t = self.threshhold
        with open(self.file_path, 'r') as f:
            samples = list(next(f).rstrip().split(self.delim)[1:])
            self.init_nearest_neighbor(samples)
            for line in f:
                line_split = line.strip().split(self.delim)
                label = line_split[0]
                self.sample_ids.append(label)

                is_found = False
                for i in range(0,len(self.components)):
                    if label in self.components[i]:
                        sample_index = i
                        is_found = True
                        break
                if not is_found:
                    self.components.append(set())
                    sample_index = -1

                distances = [float(x) for x in line_split[1:]]
                for idx,value in enumerate(distances):
                    sample_id = samples[idx]
                    if value <= t:
                        self.components[sample_index].add(sample_id)
                        d = self.near_neighbors[label]['dist']
                        if label != sample_id:
                            if d is None or d > value:
                                self.near_neighbors[label]['dist'] =value
                                self.near_neighbors[label]['sample_ids'] = set(sample_id)
                            elif d == value:
                                self.near_neighbors[label]['sample_ids'].add(sample_id)

                            d = self.near_neighbors[sample_id]['dist']
                            if d is None or d > value:
                                self.near_neighbors[sample_id]['dist'] =value
                                self.near_neighbors[sample_id]['sample_ids'] = set(label)
                            elif d == value:
                                self.near_neighbors[sample_id]['sample_ids'].add(label)


    def parse_ltriangle_distance_matrix(self):
        t = self.threshhold
        with open(self.file_path, 'r') as f:
            for line in f:
                line_split = line.strip().split(self.delim)
                label = line_split[0]
                line_split = line_split[1:]
                self.sample_ids.append(label)
                self.near_neighbors[label] = {
                            'dist':None,
                            'sample_ids':set()
                }
                is_found = False
                for i in range(0,len(self.components)):
                    if label in self.components[i]:
                        sample_index = i
                        is_found = True
                        break
                if not is_found:
                    self.components.append(set())
                    sample_index = -1
                for idx,value in enumerate(line_split):
                    if value == '':
                        break
                    sample_id = self.sample_ids[idx]

                    if float(value) <= t:
                        self.components[sample_index].add(sample_id)
                        d = self.near_neighbors[label]['dist']
                        if label != sample_id:
                            if d is None or d > value:
                                self.near_neighbors[label]['dist'] = value
                                self.near_neighbors[label]['sample_ids'] = set(sample_id)
                            elif d == value:
                                self.near_neighbors[label]['sample_ids'].add(sample_id)

                            d = self.near_neighbors[sample_id]['dist']
                            if d is None or d > value:
                                self.near_neighbors[sample_id]['dist'] = value
                                self.near_neighbors[sample_id]['sample_ids'] = set(label)
                            elif d == value:
                                self.near_neighbors[sample_id]['sample_ids'].add(label)

    def link(self):
        connected = []
        assigned = set()
        for i in range(0,len(self.components)):
            c1 = self.components[i]
            if len(c1 - assigned) == 0:
                continue
            for k in range(i+1,len(self.components)):
                c2 = self.components[k]
                ovl = c1 & c2
                if len(ovl) > 0:
                    c1 = c1 | c2
                assigned = assigned | c1
            connected.append(c1)
        self.components = connected

    def write_links(self):
        fh = open(self.out_path,"w")
        fh.write("component_id\tsample_ids\n")
        for i in range(0,len(self.components)):
            fh.write("{}\t{}\n".format(i,",".join(sorted([str(x) for x in self.components[i]]))))
        fh.close()



#n = neighbours('snp.dists.txt','snp.dists.slink.100.txt',is_square=True,thresh=100,delim='\t')
#print(n.get_result())