import pysam

def get_bamfile(bam):
    bmode = 'rb'
    if bam.endswith(('.sam', '.SAM')):
        bmode = 'r'
    elif bam.endswith(('.cram', '.CRAM')):
        bmode = 'rc'
    return pysam.AlignmentFile(bam, bmode)

def get_output(output, template):
    if output is not None:
        wbmode = 'wb'
        if output.endswith(('.sam', '.SAM')):
            wbmode = 'w'
        elif output.endswith(('.cram', '.CRAM')):
            wbmode = 'wc'
        return pysam.AlignmentFile(output, wbmode, template=template)
    return  pysam.AlignmentFile('-', "w", template=template)

