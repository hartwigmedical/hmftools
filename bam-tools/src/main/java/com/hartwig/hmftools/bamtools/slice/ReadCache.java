package com.hartwig.hmftools.bamtools.slice;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;

public class ReadCache
{
    private static final int REGION_CONSOLIDATION_DISTANCE = 10_000;

    private final Map<String,Fragment> mFragmentMap;
    private final SliceWriter mSliceWriter;
    private volatile boolean mProcessingRemoteRegions;

    public ReadCache(final SliceWriter sliceWriter)
    {
        mSliceWriter = sliceWriter;

        mFragmentMap = Maps.newHashMap();
        mProcessingRemoteRegions = false;
    }

    public void setProcessingRemoteRegions(boolean processingMateRegions)
    {
        mProcessingRemoteRegions = processingMateRegions;
    }

    public synchronized boolean addReadRecord(SAMRecord read)
    {
        if(read.isSecondaryAlignment())
        {
            if(mProcessingRemoteRegions)
                return false;

            mSliceWriter.writeRead(read); // write without tracking of fragment info
            return false;
        }

        Fragment fragment = mFragmentMap.get(read.getReadName());

        boolean isNewRead = false;

        if(fragment == null)
        {
            if(mProcessingRemoteRegions)
            {
                // don't handle new fragments after initial slice
                return false;
            }

            fragment = new Fragment(read);

            if(!fragment.isComplete())
                mFragmentMap.put(read.getReadName(), fragment);

            isNewRead = true;
        }
        else
        {
            isNewRead = fragment.processRead(read);

            if(fragment.isComplete())
                mFragmentMap.remove(read.getReadName());
        }

        if(isNewRead)
        {
            mSliceWriter.writeRead(read);
        }

        return mFragmentMap.isEmpty();
    }

    public synchronized boolean allComplete() { return mFragmentMap.isEmpty(); }

    public List<ChrBaseRegion> collateRemoteReadRegions()
    {
        Map<String,List<BaseRegion>> chrBaseRegions = Maps.newHashMap();

        for(Fragment fragment : mFragmentMap.values())
        {
            List<BasePosition> remotePositions = fragment.extractPendingRegions();

            for(BasePosition region : remotePositions)
            {
                List<BaseRegion> baseRegions = chrBaseRegions.get(region.Chromosome);

                if(baseRegions == null)
                {
                    baseRegions = Lists.newArrayList();
                    chrBaseRegions.put(region.Chromosome, baseRegions);
                }

                baseRegions.add(new BaseRegion(region.Position, region.Position));
            }
        }

        List<ChrBaseRegion> remoteRegions = Lists.newArrayList();

        // merge regions within proximity of each other
        for(Map.Entry<String,List<BaseRegion>> entry : chrBaseRegions.entrySet())
        {
            String chromosome = entry.getKey();

            List<BaseRegion> baseRegions = entry.getValue();

            Collections.sort(baseRegions);

            int index = 0;
            while(index < baseRegions.size() - 1)
            {
                BaseRegion region = baseRegions.get(index);

                int nextIndex = index + 1;
                while(nextIndex < baseRegions.size())
                {
                    BaseRegion nextRegion = baseRegions.get(nextIndex);

                    if(nextRegion.start() > region.end() + REGION_CONSOLIDATION_DISTANCE)
                        break;

                    if(nextRegion.start() <= region.end() + REGION_CONSOLIDATION_DISTANCE)
                    {
                        region.setEnd(nextRegion.end());
                        baseRegions.remove(nextIndex);
                    }
                    else
                    {
                        ++nextIndex;
                    }
                }

                ++index;
            }

            for(BaseRegion region : baseRegions)
            {
                remoteRegions.add(new ChrBaseRegion(chromosome, region.start(), region.end()));
            }
        }

        return remoteRegions;
    }

    public void logMissingReads(List<ChrBaseRegion> excludedRegions)
    {
        for(Fragment fragment : mFragmentMap.values())
        {
            List<ReadInfo> pendingReads = fragment.pendingReads();

            if(pendingReads.isEmpty())
                continue;

            boolean inExcluded = false;
            for(ChrBaseRegion excludedRegion : excludedRegions)
            {
                if(pendingReads.stream().anyMatch(x -> excludedRegion.overlaps(x.Contig, x.AlignmentStart, x.AlignmentStart)))
                {
                    inExcluded = true;
                    break;
                }
            }

            if(inExcluded)
                continue;

            BT_LOGGER.warn("read id({}) has {} missing reads:", fragment.readId(), pendingReads.size());

            for(ReadInfo readInfo : fragment.pendingReads())
            {
                BT_LOGGER.warn("missing read id({}) info({})", fragment.readId(), readInfo);
            }

            for(ReadInfo readInfo : fragment.receivedReads())
            {
                BT_LOGGER.debug("received read id({}) info({})", fragment.readId(), readInfo);
            }
        }
    }

    public String toString() { return format("fragments(%d)", mFragmentMap.size()); }

    @VisibleForTesting
    public Map<String,Fragment> fragmentMap() { return mFragmentMap; }
}
