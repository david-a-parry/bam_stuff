import pysam

def get_coordinate(read):
    if read.reference_id >= 0:
        coord = "{}:{:,}".format(read.reference_name, read.reference_start)
    else:
        coord = "*/*"
    return coord

def get_rname(read):
    if read.reference_id >= 0:
        return read.reference_name
    return '*'

