package com.hartwig.hmftools.bamtools.slice;

public enum WriteType
{
    BAM,
    READS;

    public String extension()
    {
        switch(this)
        {
            case BAM:
                return ".bam";
            case READS:
                return ".reads.tsv";
            default:
                return "unknown";
        }
    }
}
