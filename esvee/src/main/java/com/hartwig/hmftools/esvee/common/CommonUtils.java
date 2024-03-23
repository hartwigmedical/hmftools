package com.hartwig.hmftools.esvee.common;

import com.hartwig.hmftools.common.bam.BamSlicer;

public final class CommonUtils
{
    public static BamSlicer createBamSlicer()
    {
        BamSlicer bamSlicer = new BamSlicer(0, false, true, false);
        bamSlicer.setKeepUnmapped();
        bamSlicer.setKeepHardClippedSecondaries();
        return bamSlicer;
    }

}
