package com.hartwig.hmftools.svprep;

import static java.lang.Math.round;

import com.hartwig.hmftools.common.samtools.BamSlicer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvCommon
{
    public static final Logger SV_LOGGER = LogManager.getLogger(SvCommon.class);

    private static final long MEGABYTE = 1024L * 1024L;

    public static BamSlicer createBamSlicer()
    {
        BamSlicer bamSlicer = new BamSlicer(0, false, true, false);
        bamSlicer.setKeepUnmapped();
        bamSlicer.setKeepHardClippedSecondaries();
        return bamSlicer;
    }
}
