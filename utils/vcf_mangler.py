from utils.utils import parse_FORMAT
from utils import bed_reader as br


def get_newlines(line_elements, history, fasta, dp, pl):
    """
    Construct the  gVCF portions preceding this variant, up to this variant and return a new minimum value

    :param dp:
    :param pl:
    :param line_elements: variant line from VCF file
    :param history: history dictionary that maintains tracking information
    :param fasta: pytabix FASTA object
    :return: g_line, v_line, history
    """
    chrom = str(line_elements[0])
    previous_stopping_point = history['current_min']
    prev_ref_base = fasta[chrom][int(previous_stopping_point)].seq
    assert prev_ref_base != "", f"Varaint position ({chrom}:{previous_stopping_point}) not found"

    #########################################################################################
    # build g_line
    #########################################################################################
    num_samples = 1
    line = list()
    line.append(chrom)  # CHROM
    line.append(previous_stopping_point)  # POS
    line.append('.')  # ID
    line.append(prev_ref_base)  # REF
    line.append("<NON_REF>")  # ALT
    line.append('.')  # QUAL
    line.append('.')  # FILTER
    line.append('Inferred;END=' + str(int(line_elements[1]) - 1) + ';')  # INFO
    line.append("GT:MIN_DP:PL")  # FORMAT

    SAMPLE = ["0/0:" + str(dp) + ":0," + str(pl) + "," + str(pl)] * num_samples
    line.extend(SAMPLE)
    g_line = '\t'.join([str(x) for x in line])

    #########################################################################################
    # build v_line
    #########################################################################################
    new_stopping_point = line_elements[1]

    line = list()
    line.append(chrom)  # CHROM
    line.append(line_elements[1])  # POS
    line.append(line_elements[2])  # ID
    line.append(line_elements[3])  # REF
    line.append(line_elements[4] + ",<NON_REF>")  # ALT
    line.append(line_elements[5])  # QUAL
    line.append(line_elements[4])  # FILTER
    line.append('Inferred;' + line_elements[7])  # INFO
    line.append(line_elements[8])  # FORMAT

    NEW_SAMPLE = parse_FORMAT(line_elements[8], line_elements[9], pl)
    line.append(NEW_SAMPLE)
    v_line = '\t'.join([str(x) for x in line])

    #########################################################################################
    # Update history
    #########################################################################################
    history['current_max'] = int(new_stopping_point) + 1
    history['current_min'] = int(new_stopping_point)

    # Change the list of able bed regions
    history['current_regions'] = br.filter_chrom_pos_data(chrom, int(new_stopping_point) + 1, history)
    return g_line, v_line, history


def get_final_line(history, fasta, dp, pl):
    final_lines = []

    for i in history['current_regions']:
        #########################################################################################
        # build g_line
        #########################################################################################
        num_samples = 1
        prev_ref_base = fasta[i[0]][int(i[1])].seq
        line = list()
        line.append(i[0])  # CHROM
        line.append(i[1])  # POS
        line.append('.')  # ID
        line.append(prev_ref_base)  # REF
        line.append("<NON_REF>")  # ALT
        line.append('.')  # QUAL
        line.append('.')  # FILTER
        line.append('Inferred;END=' + str(i[2]) + ';')  # INFO
        line.append("GT:MIN_DP:PL")  # FORMAT

        SAMPLE = ["0/0:" + str(dp) + ":0," + str(pl) + "," + str(pl)] * num_samples
        line.extend(SAMPLE)
        g_line = '\t'.join([str(x) for x in line])
        final_lines.append(g_line)
    return final_lines
