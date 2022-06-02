package com.hartwig.hmftools.sage.common;

import java.util.List;
import java.util.function.Consumer;

import org.apache.commons.compress.utils.Lists;

import htsjdk.samtools.SAMRecord;

public class MockSamSlicer implements SamSlicerInterface
{
    public final List<SAMRecord> ReadRecords;

    public MockSamSlicer()
    {
        ReadRecords = Lists.newArrayList();
    }

    @Override
    public void slice(final Consumer<SAMRecord> consumer)
    {
        ReadRecords.forEach(x -> consumer.accept(x));
    }
}
