package com.hartwig.hmftools.sieve.annotate;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public interface IJobRegion
{
    ChrBaseRegion getChrBaseRegion();

    void matchedRead(@NotNull final SAMRecord read);
}
