package com.hartwig.hmftools.redux.common;

import java.util.List;

import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.umi.UmiConfig;
import com.hartwig.hmftools.redux.umi.UmiGroupBuilder;

import htsjdk.samtools.SAMRecord;

public class DuplicateGroupBuilder
{
    private final UmiConfig mUmiConfig;
    private final Statistics mStats;
    private final boolean mFormConsensus;
    private final UmiGroupBuilder mUmiGroupBuilder;

    public DuplicateGroupBuilder(final ReduxConfig config)
    {
        mFormConsensus = config.FormConsensus;
        mUmiConfig = config.UMIs;
        mStats = new Statistics();
        mUmiGroupBuilder = new UmiGroupBuilder(config.UMIs, mStats.UmiStats);
    }

    public Statistics statistics() { return mStats; }

    public List<DuplicateGroup> processDuplicateGroups(
            final List<DuplicateGroup> rawDuplicateGroups, final List<ReadInfo> singleFragments, boolean captureStats)
    {
        if(mUmiConfig.Enabled)
        {
            List<DuplicateGroup> umiGroups = mUmiGroupBuilder.processUmiGroups(rawDuplicateGroups, singleFragments, captureStats);

            if(captureStats)
            {
                for(DuplicateGroup umiGroup: umiGroups)
                {
                    mStats.addUmiGroup(umiGroup.readCount(), umiGroup.hasDualStrand());
                }
            }

            return umiGroups;
        }

        if(captureStats)
        {
            for(DuplicateGroup duplicateGroup : rawDuplicateGroups)
            {
                mStats.addDuplicateGroup(duplicateGroup.readCount());

                if(!mFormConsensus)
                {
                    SAMRecord primaryRead = getPrimaryRead(duplicateGroup.reads());
                    duplicateGroup.setPrimaryRead(primaryRead);
                }
            }
        }

        return rawDuplicateGroups;
    }

    public void processAllDuplicateGroups(
            final List<DuplicateGroup> duplicateGroups, final List<MultiCoordsDuplicateGroup> multiCoordsDuplicateGroups,
            boolean captureStats)
    {
        if(captureStats)
        {
            for(DuplicateGroup duplicateGroup : duplicateGroups)
            {
                mStats.addDuplicateGroup(duplicateGroup.readCount());

                if(!mFormConsensus)
                {
                    SAMRecord primaryRead = getPrimaryRead(duplicateGroup.reads());
                    duplicateGroup.setPrimaryRead(primaryRead);
                }
            }

            for(MultiCoordsDuplicateGroup multiCoordsDuplicateGroup : multiCoordsDuplicateGroups)
            {
                mStats.addDuplicateGroup(multiCoordsDuplicateGroup.readCount());

                if(!mFormConsensus)
                {
                    SAMRecord primaryRead = getPrimaryRead(multiCoordsDuplicateGroup.reads());
                    multiCoordsDuplicateGroup.setPrimaryRead(primaryRead);
                }
            }
        }
    }

    private static SAMRecord getPrimaryRead(final List<SAMRecord> reads)
    {
        SAMRecord maxRead = null;
        double maxBaseQual = 0;

        for(SAMRecord read : reads)
        {
            double avgBaseQual = calcBaseQualAverage(read);

            if(avgBaseQual > maxBaseQual)
            {
                maxBaseQual = avgBaseQual;
                maxRead = read;
            }
        }

        return maxRead;
    }

    public static double calcBaseQualAverage(final SAMRecord read)
    {
        int readBaseCount = 0;
        int readBaseQualTotal = 0;

        for(int i = 0; i < read.getBaseQualities().length; ++i)
        {
            ++readBaseCount;
            readBaseQualTotal += read.getBaseQualities()[i];
        }

        return readBaseCount > 0 ? readBaseQualTotal / (double)readBaseCount : 0;
    }
}
