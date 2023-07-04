import os

class matrix_splitter:
    file_path = None
    out_path = None
    batch_size = None
    prefix = 'segment'
    num_lines = 0
    delim = "\t"
    is_ok = True
    partitions = []
    error_msg = []

    def __init__(self,file_path,out_path,batch_size,partitions=[],prefix=None,delim='\t'):
        self.file_path = file_path
        self.out_path = out_path
        self.batch_size = batch_size
        self.prefix = prefix
        self.delim = delim
        self.num_batches = 0
        self.ranges = []
        if len(partitions) > 0:
            self.partitions = partitions

        if not os.path.isfile(file_path):
            self.error_msg.append("Error matrix file: {} does not exist".format(self.file_path))
            self.is_ok = False
            return

        self.num_lines = self.get_file_length() - 1

        if self.num_lines < 2:
            self.error_msg.append("Error matrix file: {} does not contain at least two samples".format(self.file_path))
            self.is_ok = False
            return

        if not os.path.isdir(self.out_path):
            self.error_msg.append("Directory: {} does not exist".format(self.out_path))
            self.is_ok = False
            return

        if self.batch_size  < 5:
            self.batch_size = 5


    def get_file_length(self):
        return int(os.popen(f'wc -l {self.file_path}').read().split()[0])

    def prep_batch_ranges(self):
        self.num_batches = int(self.num_lines / self.batch_size)
        rem = self.batch_size % self.num_batches
        ranges = []
        for i in range(0,self.num_batches):
            ranges.append(i*self.batch_size,i*self.batch_size+self.batch_size)
        if rem != 0:
            r = ranges[-1]
            r[1] = self.num_lines
        self.ranges = ranges


    def parse_distance_matrix_bins(self):
        '''
        Reads in a lower triangle/full distance matrix and splits it into component matricies
        according to the desired number of samples in each batch. Matrix is returned in lower triangle format
        :return:
        '''
        with open(self.file_path, 'r') as f:
            header = next(f).split(self.delim)  # skip header
            line_num = 0
            range_index = 0
            start, end = self.ranges[range_index]
            out_fh = open(os.path.join(self.out_path,"{}-{}.matrix".format(self.prefix,range_index)),'w')
            out_fh.write("{}\n".format("{}".format(self.delim).join([str(x) for x in header[start,end]])))
            for line in f:
                line_split = line.strip().split(self.delim)
                label = line_split[0]
                distances = list(map(float, line_split[start:end]))
                out_fh.write("{}\t{}\n".format(label,"{}".format(self.delim).join([str(x) for x in distances])))
                line_num += 1
                if line_num > end:
                    range_index+=1
                    start, end = self.ranges[range_index]
                    out_fh.close()
                    out_fh = open(os.path.join(self.out_path, "{}-{}.matrix".format(self.prefix, range_index)), 'w')
                    out_fh.write("{}\n".format("{}".format(self.delim).join([str(x) for x in header[start, end]])))
        out_fh.close()

    def parse_distance_matrix_partitions(self):
        '''
        Reads in a lower triangle/full distance matrix and splits it into component matricies
        according to the desired number of samples in each batch. Matrix is returned in lower triangle format
        :return:
        '''
        with open(self.file_path, 'r') as f:
            header = next(f).split(self.delim)  # skip header
            line_num = 0
            range_index = 0
            start, end = self.ranges[range_index]
            out_fh = open(os.path.join(self.out_path,"{}-{}.matrix".format(self.prefix,range_index)),'w')
            out_fh.write("{}\n".format("{}".format(self.delim).join([str(x) for x in header[start,end]])))
            for line in f:
                line_split = line.strip().split(self.delim)
                label = line_split[0]
                distances = list(map(float, line_split[start:end]))
                out_fh.write("{}\t{}\n".format(label,"{}".format(self.delim).join([str(x) for x in distances])))
                line_num += 1
                if line_num > end:
                    range_index+=1
                    start, end = self.ranges[range_index]
                    out_fh.close()
                    out_fh = open(os.path.join(self.out_path, "{}-{}.matrix".format(self.prefix, range_index)), 'w')
                    out_fh.write("{}\n".format("{}".format(self.delim).join([str(x) for x in header[start, end]])))
        out_fh.close()