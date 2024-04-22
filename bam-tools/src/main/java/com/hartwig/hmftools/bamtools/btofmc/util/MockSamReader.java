package com.hartwig.hmftools.bamtools.btofmc.util;

import java.util.Collection;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public abstract class MockSamReader
{
    public static SamReader create(final Collection<SAMRecord> samRecords)
    {
        return new SamReader.PrimitiveSamReaderToSamReaderAdapter(new MockPrimitiveSamReader(samRecords), null);
    }
}
