def parse_bed(bed_file):
    """
    Read bed file
    :param bed_file:
    :return: list of bed coord list i.e. [['20', '0', '19'], ['20', '22', '25']]
    """
    f = open(bed_file, "r")
    targets = []
    for x in f:
        d = x.strip().split('\t', 4)
        targets.append(d)
    return targets


def greater_than(value1, dict1):
    """
    Verifying that bed files are sorted by ensuring this chromosome hasn't been see before, but if it has, make sure
    end position is > the last position seen.

    :param value1: bed range [1 100 105]
    :param dict1: dictionary to keep track of what has been seen.
                        { "1": 85 } # Previously saw only up to position 85, so this would pass
                        { "1": 185 } # Previously saw up to position 185, so this would fail
    :return: dict of new max for that chrom
    """
    assert (len(value1) == 3)
    chr, start, end = value1

    # Add the chrom if I haven't seen it before
    if not str(chr) in dict1.keys():
        dict1[str(chr)] = 0

    if int(end) < int(dict1[str(chr)]):
        raise ValueError(f"Your bed file is not sorted! Previously saw up to {dict1[str(chr)]} on chrom: "
                         f"{chr} but you provided {end}")

    # Update largest event seen
    dict1[str(chr)] = end
    return dict1


def get_chrom_data(chrom, chrom_list):
    """
    Get bed regions for a particular chromosome
    :param chrom: a chromosome to extract
    :param chrom_list: list of bed coord [['20', '0', '19'], ['20', '22', '25']]
    :return: a list of regions in the chromosome of interest
    """
    return [x for x in chrom_list if x[0] == chrom]


def filter_chrom_data(chrom, chrom_list):
    """
    Exclude bed regions for a particular chromosome
    :param chrom: a chromosome to extract
    :param chrom_list: list of bed coord [['20', '0', '19'], ['20', '22', '25']]
    :return: a list of regions in the chromosome of interest
    """
    return [x for x in chrom_list if x[0] != chrom]


def filter_chrom_pos_data(chrom, pos, history):
    """
    Exclude bed regions for a particular chromosome
    :param history:
    :param chrom: a chromosome to extract
    :param pos: a minimum position
st that remain

    Ex:
        j = filter_chrom_pos_data('20', 20, [['20', '0', '19'], ['20', '22', '25']])
        ['20', '22', '25']
        j = filter_chrom_pos_data('20', 5, [['20', '0', '19'], ['20', '22', '25']])
        [['20', '5', '19'], ['20', '22', '25']]

    """
    bed_regions = []
    for x in history['current_regions']:
        if x[0] == chrom:
            # Keep if the region is within the defined region, but update the min position
            if int(x[1]) > int(pos) > int(x[2]):
                new_list = [chrom, int(pos) + 1, x[2]]
                bed_regions.append(new_list)
            # Otherwise, keep if regions are larger than the position
            elif int(x[1]) > int(pos):
                bed_regions.append(x)
    return bed_regions


def get_chrom_regions(chrom, bed_regions):
    """
    Filter the BED regions to only those for the chromosome while also returning the min, max, total of the target
    :param chrom: chromosome
    :param bed_regions: available target regions
    :return: current_max, current_min, current_regions, bed_regions
    """
    current_regions = get_chrom_data(chrom, bed_regions)
    current_max = current_regions[-1][2]
    current_min = current_regions[0][1]
    bed_regions = filter_chrom_data(chrom, bed_regions)
    return current_max, current_min, current_regions, bed_regions
