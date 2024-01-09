package com.hartwig.hmftools.sage.common;

import java.util.function.Consumer;

import htsjdk.samtools.SAMRecord;

@FunctionalInterface
public interface SamSlicerInterface
{
    void slice(final Consumer<SAMRecord> consumer);
}
