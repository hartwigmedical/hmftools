package com.hartwig.hmftools.common.bam;

public enum BamToolName
{
    SAMTOOLS,
    SAMBAMBA;

    protected static final String SAMBAMBA_COMMAND = "sambamba";
    protected static final String SAMTOOLS_COMMAND = "samtools";

    protected String toolThreadArgument()
    {
        return this == SAMTOOLS ? "-@" : "-t";
    }
}
