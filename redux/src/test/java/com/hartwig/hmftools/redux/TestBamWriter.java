package com.hartwig.hmftools.redux;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.redux.write.BamWriter;
import com.hartwig.hmftools.redux.write.ReadDataWriter;

import htsjdk.samtools.SAMRecord;

public class TestBamWriter extends BamWriter
{
    public final List<SAMRecord> WrittenRecords;
    public String CurrentChromosome;
    public int CurrentPosLower;
    public int CurrentPosUpper;

    public TestBamWriter(final ReduxConfig config)
    {
        super("", config, new ReadDataWriter(config), null, null);

        WrittenRecords = Lists.newArrayList();
        CurrentPosUpper = 0;
        CurrentPosLower = 0;
        CurrentChromosome = "";
    }

    public boolean isSorted() { return false; }

    public void initialiseRegion(final String chromosome, int startPosition)
    {
        CurrentPosLower = startPosition;
        CurrentChromosome = chromosome;
    }

    public void setBoundaryPosition(int position, boolean isLower)
    {
        if(isLower)
            CurrentPosLower = position;
        else
            CurrentPosUpper = position;
    }

    public void onRegionComplete() {}

    @Override
    protected void writeRecord(final SAMRecord read) { WrittenRecords.add(read); }

    public void close() {}

    @VisibleForTesting
    public void resetRecordWriteCounts()
    {
        mConsensusReadCount.set(0);
        mNonConsensusReadCount.set(0);
    }

}
