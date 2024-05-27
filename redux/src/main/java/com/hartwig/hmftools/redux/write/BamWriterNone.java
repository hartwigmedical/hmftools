package com.hartwig.hmftools.redux.write;

import com.hartwig.hmftools.common.basequal.jitter.JitterAnalyser;
import com.hartwig.hmftools.redux.ReduxConfig;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public class BamWriterNone extends BamWriter
{
    public BamWriterNone(final String filename, final ReduxConfig config, final ReadDataWriter readDataWriter,
            @Nullable final JitterAnalyser jitterAnalyser)
    {
        super(filename, config, readDataWriter, null, jitterAnalyser);
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
