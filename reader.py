import sys


class IO:
    def __init__(self,
                 path,
                 genome_name,
                 genome_path,
                 output_path,
                 bed_path=None):
        self.pathList = path
        self.genome_name = genome_name
        self.genome_path = genome_path
        self.bed_path = bed_path
        self.output_path = output_path

    def ref_seq(self):
        """ Returns the reference genome.
        Args:
            filename: name of the .fa reference file
        Returns:
            description: the first line of the .fa file containing the
                         description of the sequence
            short_name: short name of the sequence (needed in output file)
            result: the reference DNA sequence that has to be indexed
        """
        filename = self.genome_path
        file_object = open(filename, "r")
        ref_name = next(file_object, None)
        assert (ref_name[0] == '>'
                )  # if file in right format should start with '>'
        l = next(file_object, None)
        result = ""
        while (l != None):
            result += l[0:-1]  # don't want the '\n' symbol
            l = next(file_object, None)
        description = ref_name[1:-1]
        short_name = description.split()[0]
        # don't want the '>' and \n' symbols in the ref name
        return (description, short_name, result)

    def return_reads(self):
        reads = []
        for file in self.pathList:
            openedFile = open(file, 'r')
            for line in openedFile:
                if not (line.startswith("@")):
                    temp_line = line.split("\t")
                    if temp_line[2] == self.genome_name and temp_line[5] != '*':
                        chromosome = temp_line[0].split("-")[0]
                        reads.append([
                            chromosome, temp_line[1], temp_line[3], temp_line[5],
                            temp_line[9], temp_line[10]
                        ])

        return reads

    def read_bed(self):
        file = self.bed_path
        if file is None:
            return None
        # Returns a dict of position to variant name
        if file[-4:].lower() == '.bed':
            reads = dict()
            with open(file, 'r') as f:
                for line in f:
                    if line.startswith(self.genome_name):
                        temp_line = line.split("\t")
                        if str(int(temp_line[1]) + 1) not in reads:
                            reads[str(int(temp_line[1]) + 1)] = {}
                        reads[str(int(temp_line[1]) + 1)][int(temp_line[2]) -
                                      int(temp_line[1])] = temp_line[3]

            return reads  # returning a list to make it uniform with .sam reads

        print("Unrecognizable file ending. File endings must be .sam or .bed")
        return None

    # Columns are consisted of: (7 and 8 are not included in the output_list)
    # 0. chrom
    # 1. Pos
    # 2. ID
    # 3. reference allele
    # 4. read alleles
    # -. QUAL
    # -. FILTER
    # 5. INFO
    # 6. FORMAT
    # 7. simu
    def produce_to_file(self, output_list, chr_name, info):
        output_list = sorted(output_list, key=lambda x: x[0])
        with open(self.output_path, 'w+') as f:
            f.write("##fileformat=VCFv4.3\n")
            f.write("##reference={}\n".format(self.genome_path))
            f.write(
                '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">"\n')
            f.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" +
                    "\tN\n")
            for line in output_list:
                out_pos, name, ref_allele, base, info, gt = line
                base_list = base.split(",")
                if len(base_list) == 1:
                    if info[0][1] == "":
                        info_str = "SVLEN={}".format(info[0][0])
                    else:
                        info_str = "SVTYPE={};SVLEN={}".format(
                            info[0][1], info[0][0])
                else:
                    info_str = []
                    for i, v in enumerate(base_list):
                        info_str.append(str(info[i][0]))
                    info_str = "SVLEN=" + ",".join(info_str)

                output = "{0}\t{1}\t{2}\t{3}\t{4}\t.\t.\t{5}\t{6}\t{7}\n".format(
                    chr_name, out_pos, name, ref_allele, base, info_str, 'GT',
                    gt)
                f.write(output)

    # if __name__ == '__main__':
    #     reader = Reader('data/alignedreads.sam', 'chr22')
    #     reads = reader.return_reads()
    #     print(reads[:10])
    #     output = Output('output.vcf')
    #     output.produce_to_file([[22, 891980, 'N', 'G', 'SVLEN=1', 'GT\t0|1']])
