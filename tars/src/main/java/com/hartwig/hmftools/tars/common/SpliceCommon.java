package com.hartwig.hmftools.tars.common;

import static com.hartwig.hmftools.tars.common.TarsConstants.ALT_CONTIG_SUFFIX;

public final class SpliceCommon
{
    public static String altContigName(final String chromosome)
    {
        return chromosome + ALT_CONTIG_SUFFIX;
    }
}
