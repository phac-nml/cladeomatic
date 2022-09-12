import os,copy

class vcfReader:
    required_fields = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    row_num = 0
    row = None
    reader = None
    header = []
    samples = []

    def next_row(self):
        try:
            return next(self.reader)
        except:
            return None

    def check_file(self, file):
        if not os.path.isfile(file):
            return False
        if os.path.getsize(file) == 0:
            return False
        return True

    def __init__(self,file):

        if not self.check_file(file):
            return None

        self.reader = open(file,'r')

        while len(self.header) == 0:
            self.row = self.next_row()
            if self.row == None:
                continue
            self.row.rstrip()
            if self.row[0:6] == '#CHROM':
                self.header = self.row.rstrip().split("\t")

        for idx,value in enumerate(self.header):
            if value in self.required_fields:
                continue
            self.samples.append(str(value))

    def process_row(self):
        self.row = self.next_row()
        if self.row is None:
            return None
        line = self.row.rstrip().split("\t")
        data = {}
        for idx,value in enumerate(line):
            data[self.header[idx]] = value

        chr = data['#CHROM']
        pos = data['POS']
        ref_base = data['REF']
        alt_bases = data['ALT'].split(',')
        bases = [ref_base] + alt_bases

        for sample_id in self.samples:
            base = data[sample_id]
            if isinstance(base,str):
                if base.isnumeric():
                    base = bases[int(base)]
                elif base not in ['A', 'T', 'C', 'G','-']:
                    if base == '*':
                        base = '-'
                    else:
                        base = 'N'
            else:
                base = bases[int(base)]
            data[sample_id] = base

        return data









