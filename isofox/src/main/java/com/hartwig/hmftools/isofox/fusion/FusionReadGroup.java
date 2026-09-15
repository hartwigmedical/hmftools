package com.hartwig.hmftools.isofox.fusion;

import java.util.List;
import java.util.Map;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.common.Read;

public class FusionReadGroup
{
    public final String ReadId;

    private final List<FusionRead> mReads;
    private boolean mIsComplete;

    public FusionReadGroup(final String readId, final List<FusionRead> reads)
    {
        ReadId = readId;

        mReads = Lists.newArrayListWithCapacity(reads.size());
        mReads.addAll(reads);
        mIsComplete = readGroupComplete();
    }

    public int size() { return mReads.size(); }

    public List<FusionRead> reads() { return mReads; }
    public void addRead(final FusionRead read)
    {
        mReads.add(read);
        mIsComplete = readGroupComplete();
    }

    public boolean isComplete() { return mIsComplete; }

    public boolean hasSuppAlignment() { return hasSuppAlignment(mReads); }

    public static boolean hasSuppAlignment(final List<FusionRead> reads)
    {
        return reads.stream().anyMatch(x -> x.HasSuppAlignment);
    }

    public void merge(final FusionReadGroup other)
    {
        other.reads().forEach(x -> addRead(x));
    }

    public String toString()
    {
        return String.format("%s reads(%d) complete(%s)", ReadId, mReads.size(), mIsComplete);
    }

    public String findOtherChromosome(final String chromosome)
    {
        for(FusionRead read : mReads)
        {
            if(!read.MateChromosome.equals(chromosome))
                return read.MateChromosome;

            if(read.SuppData != null)
                return read.SuppData.Chromosome;
        }

        return null;
    }

    private boolean readGroupComplete()
    {
        int suppCount = 0;
        int nonSuppCount = 0;
        int expectedSuppCount = 0;
        int expectedNonSuppCount = 1;

        for(FusionRead read : mReads)
        {
            if(read.isReadPaired() && !read.isMateUnmapped())
            {
                expectedNonSuppCount = 2;
            }

            if(read.isSupplementaryAlignment())
            {
                ++suppCount;
            }
            else
            {
                ++nonSuppCount;

                if(read.HasSuppAlignment)
                {
                    ++expectedSuppCount;
                }
            }
        }

        return (expectedNonSuppCount == nonSuppCount) && (expectedSuppCount == suppCount);
    }

    public static void mergeChimericReadMaps(
            final Map<String,FusionReadGroup> partialGroups, final List<FusionReadGroup> completeGroups,
            final Map<String,FusionReadGroup> sourceMap)
    {
        // 1. copies complete groups from the source map into complete groups map
        // 2. checks for a partial match by combining partials and source, and if found removes from partials
        // 3. new partial groups from the source map are copied into the partials map
        // note: source map is logically const
        for(Map.Entry<String, FusionReadGroup> entry : sourceMap.entrySet())
        {
            FusionReadGroup srcReadGroup = entry.getValue();

            if(srcReadGroup.isComplete())
            {
                completeGroups.add(srcReadGroup);
            }
            else
            {
                // look for an existing incomplete group to add these reads to
                String readId = entry.getKey();
                FusionReadGroup existingReadGroup = partialGroups.get(readId);

                if(existingReadGroup == null)
                {
                    partialGroups.put(readId, srcReadGroup);
                }
                else
                {
                    existingReadGroup.merge(srcReadGroup);

                    if(existingReadGroup.isComplete())
                    {
                        partialGroups.remove(readId);
                        completeGroups.add(existingReadGroup);
                    }
                }
            }
        }
    }
}
