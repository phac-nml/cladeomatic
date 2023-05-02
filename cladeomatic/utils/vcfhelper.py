import os,copy
import cladeomatic.constants as CONSTANT

class vcfReader:
    """
    The vcfReader class is meant to be a helper class to read and
    process VCF files.
    """
    #the list of required headers for a VCF file
    required_fields = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    row_num = 0
    row = None
    reader = None
    header = []
    samples = []

    def next_row(self):
        """
        A method to return the next row of the VCF file.

        Returns
        -------
        vcfReader object
            The next row of the file
        """
        try:
            return next(self.reader)
        except:
            return None

    def check_file(self, file):
        """
        A helper method to verify that a file exists and is not empty.

        Parameters
        ----------
        file : str
            The path to the file
        Returns
        -------
        bool
            True if the file both exists and is not empty
        """
        if not os.path.isfile(file):
            return False
        if os.path.getsize(file) == CONSTANT.MIN_FILE_SIZE:
            return False
        return True

    def __init__(self,file):
        """
        The vcfReader class instantiation.  Also used to run the
        processing methods like file validation, reading the file,
        cleaning the file Strings, and creating the header and sample lists.

        Parameters
        ----------
        file : str
            The path to the VCF file
        """
        #return None if the file validation fails
        if not self.check_file(file):
            return None

        self.reader = open(file,'r')
        #create the header list
        while len(self.header) == 0:
            self.row = self.next_row()
            if self.row == None:
                continue
            self.row.rstrip()
            if self.row[0:6] == '#CHROM':
                self.header = self.row.rstrip().split("\t")
        #create the sample list
        for idx,value in enumerate(self.header):
            if value in self.required_fields:
                continue
            self.samples.append(str(value))

    def process_row(self):
        """
        A method to process the output of a row for the VCF file.

        Returns
        -------
        dict
            A dictionary of the VCF file data with the headers as keys for the sample data.  Return 'None' if there is no next row in the file.
        """
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









