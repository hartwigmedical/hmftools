package com.hartwig.hmftools.markdups.write;

import com.hartwig.hmftools.markdups.MarkDupsConfig;

import htsjdk.samtools.SAMRecord;

public class BamWriterNone extends BamWriter
{
    public BamWriterNone(
            final String filename, final MarkDupsConfig config, final ReadDataWriter readDataWriter)
    {
        super(filename, config, readDataWriter, null);
    }

    public boolean isSorted() { return false; }

    public void initialiseRegion(final String chromosome, int startPosition) {}
    public void setBoundaryPosition(int position, boolean isLower) {}
    public void onRegionComplete() {}

    @Override
    protected void writeRecord(final SAMRecord read) {}

    @Override
    public void close() {}
}
