# ########################################################################
# VCF FORMAT FIELD PARSING
# ########################################################################

def add_default_quals(sample_value, qual):
    sample_value + ',' + qual + ',' + qual + ',' + qual
    return sample_value


def add_zero(sample_value):
    return sample_value + ',0'


def parse_FORMAT(format_data, sample_data, qual):
    format_fields = format_data.split(':')
    sample_fields = sample_data.split(':')
    new_sample_fields = []

    for i in range(len(format_fields)):
        target = format_fields[i]
        sample_value = sample_fields[i]
        if target in ['AD']:
            sample_value = add_zero(sample_value)

        elif target in ['PL', 'GP', 'GL']:
            sample_value = add_default_quals(sample_value, qual)

        # If its not either of those, don't change the value
        new_sample_fields.append(sample_value)

    sample_list = ':'.join([str(x) for x in new_sample_fields])
    return sample_list
