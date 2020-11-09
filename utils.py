import re
import pickle


def CIGAR2list(CIGAR):
    temp = list(filter(None, re.split(r"(M|I|D|X)", CIGAR)))
    for i, x in enumerate(temp):
        try:
            temp[i] = int(x)
        except ValueError:
            pass
    return temp


def calculate_e(Q):
    prob = 1.0
    for char in Q:
        prob *= 1 - (10**-((ord(char) - 33) / 10.0))
    # probability that all characters are correct
    return 1 - prob  # probability that at least one is wrong


def phred_2_error(phred):
    return 10**-(phred / 10.0)


def countAllReads(genome, reads, var_thresh=0.85):
    # current data structure
    # d = {genome_pos:[0, (), (), (), ...], ...}
    # d['0'] returns a list for genome position 0
    # d['0'][0] number of non-variants
    # d['0'][1], d['0'][2], d['0'][3] are tuples(candidates)
    # Tuple structure: (reference letter(s), read letter(s), quality(as a letter))
    phred = 20
    error = phred_2_error(phred)
    reads = sorted(reads, key=lambda read: read[2])
    res = {}
    post_process_required = set()
    last_checked = int(reads[0][2]) - 2
    for k, temp in enumerate(reads):
        chromosome, flag, pos, CIGAR, read, qual = temp
        if k % 300 == 0:
            print("{} out of {} reads processed".format(k, len(reads)))
        current = AlignedRead(genome, pos, CIGAR, read, qual)
        temp = current.process_read()
        for gPos, value in temp.items():
            v1, v2, q = value
            if v1 == v2:
                if gPos not in res:
                    res[gPos] = [0]
                res[gPos][0] += 1
            elif calculate_e(q) < error:
                if gPos not in res:
                    res[gPos] = [0]
                res[gPos].append((v1, v2, q, k))
                if len(v1) > 1:
                    post_process_required.add(int(gPos))
        if last_checked < current.pos:
            for j in range(last_checked, current.pos):
                try:
                    perc = res[str(j)][0] / (
                        len(res[str(j)]) - 1 + res[str(j)][0])
                    if perc >= var_thresh:
                        res.pop(str(j))
                except KeyError:
                    pass
            last_checked = current.pos
    for j in range(last_checked, last_checked + 200):
        try:
            perc = res[str(j)][0] / (len(res[str(j)]) - 1 + res[str(j)][0])
            if perc >= var_thresh:
                res.pop(str(j))
        except KeyError:
            pass
    print("{} out of {} reads processed".format(len(reads), len(reads)))
    # Post processing to merge candidate variants
    post_process_ordered = sorted(list(post_process_required))
    for p in post_process_ordered:
        print("Post processing genome position {}...".format(p))
        postProcessPosition(res, reads, genome, p, post_process_required)
    return res


def find_longest_window(res, genome, pos, post_process_required):
    # assumption: pos is in post_process_required
    str_pos = str(pos)
    int_pos = int(pos)
    tuples = res[str_pos][1:]
    longest_in_ref = len(max(tuples, key=lambda x: len(x[0]))[0])
    reached = int_pos + longest_in_ref  # excluding reached
    current = int_pos + 1
    while current in post_process_required:
        # calculate max reach
        try:
            str_current = str(current)
            tuples = res[str_current][1:]
            longest_in_ref = len(max(tuples, key=lambda x: len(x[0]))[0])
            current_reached = current + longest_in_ref
            reached = max(reached, current_reached)
        except KeyError:
            pass
        current += 1
    for i in range(current, reached):
        try:
            tuples = res[str(i)][1:]
            longest_in_ref = len(max(tuples, key=lambda x: len(x[0]))[0])
            current_reached = i + longest_in_ref
            reached = max(reached, current_reached)
        except KeyError:
            pass

    return reached - int_pos


def postProcessPosition(res, reads, genome, pos, post_process_required):
    if pos in post_process_required and str(pos) not in res:
        return

    # if one variant is long enough that it reaches another variant
    # nearby, then those two entries are merged. In addition,
    # the longest reference allele is used in the REF entry in the VCF e.g:
    # At position 5: TA becomes a T
    # At position 5: T becomes a C
    # At position 6: A becomes a C
    # TA reaches into position 6.

    # Solution 1:
    # # Have two entries in the VCF file for position 5
    # # One entry has REF=TA and ALT=T,
    # # The other entry has REF=T and ALT=C
    # # Problem with this solution is that the genotype becomes ambiguous since
    # # we have split the entry to two.

    # Solution 2:
    # # With a single entry, set REF=TA and ALT=T,CA
    # # Problem with this is that we are assuming that the second base in ref
    # # remained unchanged in our entry for position 5. The change at position 6
    # # is instead stored in a second entry. Therefore we cannot capture the case
    # # where ALT=T, CC(i.e both T and A change always together)

    # Solution 3(Implemented):
    # Detect when a reference allele reaches another position and reprocess
    # the ONLY relevant reads from scratch taking this into account.

    # The final entry in the vcf is that at position 5 we have:
    # REF=TA ALT=T, CC
    win_length = find_longest_window(res, genome, pos, post_process_required)
    relevant_reads = set()
    ref_allele = genome[pos:pos + win_length]
    min_common = res[str(pos)][0]
    for current in range(pos, pos + win_length):
        # NOTE: This is a heuristic to find how many reads match the reference
        # at all positions in our window
        try:
            min_common = min(min_common, res[str(current)][0])
            for _, _, _, k in res[str(current)][1:]:
                relevant_reads.add(k)
            res.pop(str(current))
        except KeyError:
            pass

    res[str(pos)] = [min_common]
    for k in relevant_reads:
        r_chromosome, r_flag, r_pos, r_CIGAR, r_read, r_qual = reads[k]
        r = AlignedRead(genome, r_pos, r_CIGAR, r_read, r_qual, cache=True)
        acc_reads = []
        acc_quals = []
        for current in range(pos, pos + win_length):
            try:
                (_, temp_read, temp_qual) = r.readByGenomeIndex(current)
                acc_reads.append(temp_read)
                acc_quals.append(temp_qual)
            except KeyError:
                pass
        acc_reads = "".join(acc_reads)
        acc_quals = "".join(acc_quals)
        res[str(pos)].append((ref_allele, acc_reads, acc_quals, k))

    return


class AlignedRead:

    genome = None

    def __init__(self, genome, pos, CIGAR, read, qual, cache=False):
        if AlignedRead.genome is None:
            AlignedRead.genome = genome
        self.pos = int(pos) - 1
        self.CIGAR = CIGAR2list(CIGAR)
        self.read = read
        self.qual = qual
        self.cache = cache
        self.result = None

    def process_read(self):
        # TODO: Test assumption that there is never an
        # insertion following a deletion.
        if self.result is not None:
            return self.result
        genome = AlignedRead.genome
        result = {}
        g_index = self.pos
        r_index = 0
        for i in range(0, len(self.CIGAR), 2):
            num = self.CIGAR[i]
            letter = self.CIGAR[i + 1]
            if letter in ('M', 'X'):
                for j in range(num - 1):
                    result[str(g_index)] = (genome[g_index],
                                            self.read[r_index],
                                            self.qual[r_index])
                    g_index += 1
                    r_index += 1
                # process last one here looking ahead for insertions/deletions
                if i != (len(self.CIGAR) - 2):
                    nextNum = self.CIGAR[i + 2]
                    nextLetter = self.CIGAR[i + 3]
                    if nextLetter == 'I':
                        result[str(g_index)] = (
                            genome[g_index],
                            self.read[r_index:r_index + 1 + nextNum],
                            self.qual[r_index:r_index + 1 + nextNum])
                        g_index += 1
                        r_index += 1 + nextNum
                    elif nextLetter == 'D':
                        result[str(g_index)] = (
                            genome[g_index:g_index + 1 + nextNum],
                            self.read[r_index], self.qual[r_index])
                        g_index += 1 + nextNum
                        r_index += 1
                    elif nextLetter == 'X':
                        pass
                    else:
                        raise Exception("Unsupported CIGAR format: {}".format(
                            self.CIGAR))
                else:
                    result[str(g_index)] = (genome[g_index],
                                            self.read[r_index],
                                            self.qual[r_index])
                    g_index += 1
                    r_index += 1
            elif letter in ('I', 'D'):
                pass
            else:
                raise Exception("Unsupported CIGAR format: {}".format(
                    self.CIGAR))
        if self.cache is True:
            self.result = result
        return result

    def readByGenomeIndex(self, index):
        # index should be 0-based
        result = self.process_read()
        return result[str(index)]


def pickle_result(item, filename):
    afile = open(filename, 'wb')
    pickle.dump(item, afile)
    afile.close()


def open_jar(filename):
    file2 = open(filename, 'rb')
    new_d = pickle.load(file2)
    file2.close()
    return new_d
