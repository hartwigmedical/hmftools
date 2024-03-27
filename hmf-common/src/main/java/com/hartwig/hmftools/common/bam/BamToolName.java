package com.hartwig.hmftools.common.bam;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public enum BamToolName
{
    SAMTOOLS,
    SAMBAMBA;

    protected static final String SAMBAMBA_BINARY = "sambamba";
    protected static final String SAMTOOLS_BINARY = "samtools";

    public static final String SAMTOOLS_PATH = "samtools";
    public static final String SAMBAMBA_PATH = "sambamba";

    public static final String BAMTOOL_PATH = "bamtool";

    protected String toolThreadArgument()
    {
        return this == SAMTOOLS ? "-@" : "-t";
    }

    public static BamToolName fromPath(final String path)
    {
        if(path.contains(SAMBAMBA_BINARY))
            return SAMBAMBA;
        if(path.contains(SAMTOOLS_BINARY))
            return SAMTOOLS;
        else
            return null;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(SAMTOOLS_PATH, false, "Path to samtools for sort/merge/index");
        configBuilder.addPath(SAMBAMBA_PATH, false, "Path to sambamba for sort/merge");
        configBuilder.addPath(BAMTOOL_PATH, false, "Path to sambamba or samtools for sort/merge");
    }
}
