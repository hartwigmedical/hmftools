package com.hartwig.hmftools.bamtools.depth;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.depth.FinderConfig.DEFAULT_HIGH_DEPTH_REGION_MAX_GAP;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._Y;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.HighDepthRegion;
import com.hartwig.hmftools.common.utils.Integers;

public class CombinerMergeTask implements Callable<Void>
{
    private final CombinerConfig mConfig;
    private final String mChromosome;
    private List<List<HighDepthRegion>> mSampleRegions;
    private final List<CombinedRegion> mCombinedRegions;
    private final List<HighDepthRegion> mHighDepthRegions;

    public CombinerMergeTask(final CombinerConfig config, final String chromosome, final List<List<HighDepthRegion>> sampleRegions)
    {
        mConfig = config;
        mChromosome = chromosome;
        mSampleRegions = sampleRegions;
        mCombinedRegions = Lists.newArrayList();
        mHighDepthRegions = Lists.newArrayList();
    }

    public String chromosome() { return mChromosome; }
    public List<HighDepthRegion> highDepthRegions() { return mHighDepthRegions; }

    public Void call()
    {
        BT_LOGGER.info("merging chromosome({})", mChromosome);

        // first generate regions from across all samples
        mergeSampleRegions();

        // then merge regions from all samples
        mHighDepthRegions.addAll(mergeChromosomeRegions(mChromosome, mCombinedRegions));

        mCombinedRegions.clear(); // to free up memory

        BT_LOGGER.debug("chromosome({}) merge complete", mChromosome);

        return null;
    }

    private void mergeSampleRegions()
    {
        BT_LOGGER.debug("chromosome({}) merging {} sample regions", mChromosome, mSampleRegions.size());

        int sampleIndex = 1;

        for(List<HighDepthRegion> regions : mSampleRegions)
        {
            BT_LOGGER.trace("merging sample({})", sampleIndex++);

            for(HighDepthRegion region : regions)
            {
                int index = 0;
                boolean matched = false;

                while(index < mCombinedRegions.size())
                {
                    CombinedRegion combinedRegion = mCombinedRegions.get(index);
                    if(positionsOverlap(region.start(), region.end(), combinedRegion.start(), combinedRegion.end()))
                    {
                        matched = true;
                        combinedRegion.addBases(region);
                        break;
                    }
                    else if(region.end() < combinedRegion.start())
                    {
                        break;
                    }

                    ++index;
                }

                if(!matched)
                {
                    CombinedRegion combinedRegion = new CombinedRegion(region);
                    mCombinedRegions.add(index, combinedRegion);
                }
                else
                {
                    // check if this matched region now overlaps with following ones
                    CombinedRegion matchedRegion = mCombinedRegions.get(index);

                    int nextIndex = index + 1;
                    while(nextIndex < mCombinedRegions.size())
                    {
                        CombinedRegion combinedRegion = mCombinedRegions.get(nextIndex);

                        if(!positionsOverlap(matchedRegion.start(), matchedRegion.end(), combinedRegion.start(), combinedRegion.end()))
                            break;

                        matchedRegion.addRegion(combinedRegion);
                        mCombinedRegions.remove(nextIndex);
                    }
                }
            }
        }
    }

    private List<HighDepthRegion> mergeChromosomeRegions(final String chromosome, final List<CombinedRegion> combinedRegions)
    {
        BT_LOGGER.debug("chromosome({}) merging {} combined regions", mChromosome, combinedRegions.size());

        // form consolidated high depth regions from the combined regions position data
        List<HighDepthRegion> highDepthRegions = Lists.newArrayList();

        int minSampleCount = mConfig.MinSampleCount;

        if(HumanChromosome.fromString(chromosome) == _Y)
            minSampleCount *= HighDepthCombiner.CHROMOSOME_Y_SAMPLE_FRACTION;

        for(CombinedRegion region : combinedRegions)
        {
            HighDepthRegion currentRegion = null;
            List<Integer> positionAverages = Lists.newArrayList();

            for(int i = 0; i < region.Depth.size(); ++i)
            {
                PositionCount positionCount = region.Depth.get(i);

                if(positionCount.count() >= minSampleCount)
                {
                    if(currentRegion == null)
                    {
                        currentRegion = new HighDepthRegion(new ChrBaseRegion(mChromosome, positionCount.Position, positionCount.Position));
                        currentRegion.DepthMin = positionCount.DepthMin;
                        currentRegion.DepthMax = positionCount.DepthMax;
                        currentRegion.SampleCount = positionCount.count();

                        positionAverages.add(positionCount.medianDepth());

                        highDepthRegions.add(currentRegion);
                    }
                    else
                    {
                        // extend the region
                        currentRegion.setEnd(positionCount.Position);
                        currentRegion.DepthMin = min(currentRegion.DepthMin, positionCount.DepthMin);
                        currentRegion.DepthMax = max(currentRegion.DepthMax, positionCount.DepthMax);

                        currentRegion.SampleCount = max(currentRegion.SampleCount, positionCount.count());

                        positionAverages.add(positionCount.medianDepth());
                    }
                }
                else
                {
                    if(currentRegion == null)
                        continue;

                    if(positionCount.Position - currentRegion.end() < DEFAULT_HIGH_DEPTH_REGION_MAX_GAP)
                        continue;

                    // end this region
                    currentRegion.DepthAvg = (int)round(Integers.median(positionAverages));
                    positionAverages.clear();
                    currentRegion = null;
                }
            }

            if(currentRegion != null)
                currentRegion.DepthAvg = (int)round(Integers.median(positionAverages));
        }

        // check min width for the region
        if(mConfig.MinRegionSize > 0)
        {
            for(HighDepthRegion highDepthRegion : highDepthRegions)
            {
                if(highDepthRegion.baseLength() < mConfig.MinRegionSize)
                {
                    int diff = mConfig.MinRegionSize - highDepthRegion.baseLength();
                    int halfExtension = diff / 2;

                    highDepthRegion.setStart(highDepthRegion.start() - halfExtension);
                    highDepthRegion.setEnd(highDepthRegion.end() + halfExtension);
                }
            }

            ChrBaseRegion.checkMergeOverlaps(highDepthRegions, true);
        }

        return highDepthRegions;
    }
}
