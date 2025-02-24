package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.commons.lang3.Validate;

import htsjdk.samtools.SAMRecord;

public class ReadCache
{
    private static final int REGION_CONSOLIDATION_DISTANCE = 10_000;

    private final Map<String, FragmentReadTracker> mFragmentReadTrackers = new HashMap<>();
    private final SliceWriter mSliceWriter;
    private volatile boolean mProcessingRemoteRegions = false;

    public ReadCache(final SliceWriter sliceWriter)
    {
        mSliceWriter = sliceWriter;
    }

    public void setProcessingRemoteRegions(boolean processingMateRegions)
    {
        mProcessingRemoteRegions = processingMateRegions;
    }

    public synchronized void addReadRecord(SAMRecord read)
    {
        FragmentReadTracker fragmentReadTracker = mFragmentReadTrackers.get(read.getReadName());

        if(fragmentReadTracker == null)
        {
            if(mProcessingRemoteRegions)
            {
                // we do not accept new fragments after initial slice
                return;
            }

            fragmentReadTracker = new FragmentReadTracker(read.getReadName(), read.getReadPairedFlag());
            mFragmentReadTrackers.put(read.getReadName(), fragmentReadTracker);
        }
        if(fragmentReadTracker.processRead(read))
        {
            // this read is new, write it out
            mSliceWriter.writeRead(read);
        }
    }

    public List<ChrBaseRegion> getRemoteReadRegions()
    {
        List<BasePosition> remoteReadPositions = new ArrayList<>();
        for(FragmentReadTracker fragmentReadTracker : mFragmentReadTrackers.values())
        {
            assert fragmentReadTracker.invariant();
            remoteReadPositions.addAll(fragmentReadTracker.findPendingReadPositions());
        }

        // sort the mate regions
        remoteReadPositions.sort(Comparator.comparing((BasePosition o) -> o.Chromosome).thenComparingInt(o -> o.Position));
        List<ChrBaseRegion> remoteBaseRegions = new ArrayList<>();

        // now we go through the remote positions and create a new condense list
        ChrBaseRegion overlapRegion = null;
        for(BasePosition br : remoteReadPositions)
        {
            // see if still overlap
            if(overlapRegion != null && br.Chromosome.equals(overlapRegion.chromosome()))
            {
                // we should be sorted this way
                Validate.isTrue(br.Position >= overlapRegion.start());
                if(overlapRegion.end() + REGION_CONSOLIDATION_DISTANCE >= br.Position)
                {
                    // still overlapping, we can set end
                    if(br.Position > overlapRegion.end())
                    {
                        overlapRegion.setEnd(br.Position);
                    }
                    continue;
                }
            }
            // no longer can reuse old one, make a new one
            overlapRegion = new ChrBaseRegion(br.Chromosome, br.Position, br.Position);
            remoteBaseRegions.add(overlapRegion);
        }
        return remoteBaseRegions;
    }

    public void logMissingReads()
    {
        for(FragmentReadTracker fragmentReadTracker : mFragmentReadTrackers.values())
        {
            for(FragmentReadTracker.ReadKey readKey : fragmentReadTracker.getPendingReads())
            {
                BT_LOGGER.warn("read id({}) missing read({})", fragmentReadTracker.getName(), readKey);
            }
        }
    }
}
