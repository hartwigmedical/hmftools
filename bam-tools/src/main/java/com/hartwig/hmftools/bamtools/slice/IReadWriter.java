package com.hartwig.hmftools.bamtools.slice;

import htsjdk.samtools.SAMRecord;

public interface IReadWriter
{
    void writeRead(final SAMRecord record);
}
