package com.hartwig.hmftools.bamtools.slice;

public enum WriteType
{
    BAM,
    READS;

    public String extension()
    {
        return switch(this)
        {
            case BAM -> ".bam";
            case READS -> ".reads.tsv";
            default -> "unknown";
        };
    }
}
