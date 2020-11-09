import collections

import numpy as np

# Note: given input candidate positions are 0-based
# only look at variant if it's probability of wrongly mapped is less than 0.01 (phred = 20)
# Returns: a list of tuples
#           (pos, name, reference_allele, variant_alleles, simu notation e.g. 0|1)


def naive_caller(candidates, genome, bed):
    output = []
    for pos in candidates:  # pos is str key
        t_variants = candidates[pos][1:]
        if len(t_variants) > 0:
            ref_allele = t_variants[0][0]
        else:
            ref_allele = genome[int(pos)]
        variants = [base for r, base, q, _ in t_variants]
        ref_count = candidates[pos][0]
        total_count = ref_count + len(variants)

        # check for variants when ref% is < 60%
        if (len(variants) != 0) and (ref_count / total_count < 0.60):
            alleles = dict()
            for base in variants:
                if not base in alleles:
                    alleles[base] = 1
                else:
                    alleles[base] += 1

            homozygous = False
            o_var = []
            for key in alleles:
                b_percent = alleles[key] / total_count
                if b_percent > 0.7:
                    homozygous = True
                    o_var.append(key)
                    break
                elif b_percent > 0.4:
                    homozygous = False
                    o_var.append(key)

            # After collecting the most probable variants, sometimes we
            # end up in a case where (e.g) the ref=TATGCACAGA and ALT=TAGG
            # This should be modified to REF=ATGCACAGA and ALT=AGG
            offset = 0
            final_pos = pos
            if len(ref_allele) > 1 and len(o_var) > 0:
                min_variant = len(min(o_var, key=lambda x: len(x)))
                i = 0
                for i in range(0, min_variant - 1):
                    flag = True
                    for v in o_var:
                        if ref_allele[i] != v[i] or ref_allele[i + 1] != v[i +
                                                                           1]:
                            flag = False
                    if flag is False:
                        break
                offset = i

                ref_allele = ref_allele[offset:]
                final_pos = str(int(pos) + offset)
                for i, v in enumerate(o_var):
                    o_var[i] = v[offset:]

            # if o_var is empty that means ref % is low due to different read errors
            if len(o_var) != 0:
                key_list = sorted(o_var)
                key_dict = {ref_allele: 0}
                for i in range(len(key_list)):
                    key_dict[key_list[i]] = i + 1

                info = []
                for v in key_list:
                    svlen = len(v) - len(ref_allele)
                    svtype = ""
                    if svlen > 0 and v.startswith(ref_allele):
                        svtype = "INS"
                    if svlen < 0 and ref_allele.startswith(v):
                        svtype = "DEL"
                    if svlen == 0:
                        svlen = 1
                    info.append((svlen, svtype))

                out_pos = str(int(final_pos) + 1)  # output position is 1-based
                name = "NA"
                if bed is not None:
                    name = bed.get(out_pos, "NA")
                    if name == "NA":
                        continue  # if BED is given then include only those in BED
                    else:
                        index = min(
                            name, key=lambda x: abs(x - len(ref_allele)))
                        name = name[index]
                if homozygous:
                    base = key_list[0]
                    output.append((out_pos, name, ref_allele, base,
                                   info, '{}|{}'.format(
                                       key_dict[base], key_dict[base])))
                else:
                    # heterozygous but not matching ref.
                    if len(key_list) == 2:
                        base = ",".join(key_list)
                        output.append((out_pos, name, ref_allele, base, info,
                                       '{}|{}'.format(key_dict[key_list[0]],
                                                      key_dict[key_list[1]])))
                    # heterozygous but matching ref.
                    if len(key_list) == 1:
                        base = key_list[0]

                        output.append((out_pos, name, ref_allele, base, info,
                                       '{}|{}'.format(key_dict[ref_allele],
                                                      key_dict[base])))
    return output
