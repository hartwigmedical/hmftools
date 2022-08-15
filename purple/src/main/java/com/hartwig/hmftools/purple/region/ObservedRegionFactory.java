package com.hartwig.hmftools.purple.region;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.genome.window.Window;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.purple.segment.PurpleSegment;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;

public class ObservedRegionFactory
{
    private final int mWindowSize;
    private final CobaltChromosomes mCobaltChromosomes;
    private final GermlineStatusFactory mStatusFactory;

    public ObservedRegionFactory(final int windowSize, final CobaltChromosomes cobaltChromosomes)
    {
        mWindowSize = windowSize;
        mCobaltChromosomes = cobaltChromosomes;
        mStatusFactory = new GermlineStatusFactory(cobaltChromosomes);
    }

    public List<ObservedRegion> combine(
            final List<PurpleSegment> regions, final Multimap<Chromosome, AmberBAF> bafs,
            final Map<Chromosome,List<CobaltRatio>> ratios, final Multimap<Chromosome, GCProfile> gcProfiles)
    {
        final List<ObservedRegion> result = Lists.newArrayList();

        final GenomePositionSelector<CobaltRatio> cobaltSelector = GenomePositionSelectorFactory.create(ratios);
        final GenomePositionSelector<AmberBAF> bafSelector = GenomePositionSelectorFactory.create(bafs);
        final GenomeRegionSelector<GCProfile> gcSelector = GenomeRegionSelectorFactory.createImproved(gcProfiles);

        for(final PurpleSegment region : regions)
        {
            final BAFAccumulator baf = new BAFAccumulator();
            final CobaltAccumulator cobalt = new CobaltAccumulator(mWindowSize, region);
            final GCAccumulator gc = new GCAccumulator(region);

            bafSelector.select(region, baf);
            cobaltSelector.select(region, cobalt);
            gcSelector.select(region, gc);

            double tumorRatio = cobalt.tumorMedianRatio();
            // double tumorRatio = cobalt.tumorMeanRatio();
            double normalRatio = cobalt.referenceMeanRatio();
            int depthWindowCount = cobalt.tumorCount();
            GermlineStatus germlineStatus = mStatusFactory.calcStatus(region.chromosome(), normalRatio, tumorRatio, depthWindowCount);

            final ObservedRegion observedRegion = new ObservedRegion(
                    region.chromosome(), region.start(), region.end(), region.RatioSupport, region.Support, baf.count(), baf.medianBaf(),
                    depthWindowCount, tumorRatio, normalRatio, cobalt.unnormalisedReferenceMeanRatio(), germlineStatus,
                    region.SvCluster, gc.averageGCContent(), region.MinStart, region.MaxStart);

            result.add(observedRegion);
        }

        return extendMinSupport(result);
    }

    public static List<ObservedRegion> extendMinSupport(final List<ObservedRegion> modifiables)
    {
        for(int i = 0; i < modifiables.size(); i++)
        {
            final ObservedRegion target = modifiables.get(i);
            if(target.support() == SegmentSupport.NONE && target.germlineStatus() == GermlineStatus.DIPLOID)
            {
                for(int j = i - 1; j >= 0; j--)
                {
                    final ObservedRegion prior = modifiables.get(j);
                    if(prior.germlineStatus() == GermlineStatus.DIPLOID)
                    {
                        break;
                    }

                    target.setMinStart(Math.min(target.minStart(), prior.start()));

                    if(prior.support() != SegmentSupport.NONE)
                    {
                        break;
                    }
                }
            }
        }
        return new ArrayList<>(modifiables);
    }

    private class BAFAccumulator implements Consumer<AmberBAF>
    {
        private int mCount;
        final private List<Double> mBafs = Lists.newArrayList();

        @Override
        public void accept(final AmberBAF baf)
        {
            if(mCobaltChromosomes.contains(baf.chromosome()))
            {
                CobaltChromosome cobaltChromosome = mCobaltChromosomes.get(baf.chromosome());
                if(cobaltChromosome.isNormal() && cobaltChromosome.isDiploid() && !Double.isNaN(baf.tumorModifiedBAF()))
                {
                    mCount++;
                    mBafs.add(baf.tumorModifiedBAF());
                }
            }
        }

        private int count()
        {
            return mCount;
        }

        private double medianBaf()
        {
            if(mCount > 0)
            {
                Collections.sort(mBafs);
                return mBafs.size() % 2 == 0 ? (mBafs.get(mCount / 2) + mBafs.get(mCount / 2 - 1)) / 2 : mBafs.get(mCount / 2);
            }
            return 0;
        }
    }

    @VisibleForTesting
    static class CobaltAccumulator implements Consumer<CobaltRatio>
    {
        private final Window mWindow;
        private final GenomeRegion mRegion;

        private final RatioAccumulator mReferenceAccumulator;
        private final RatioAccumulator mUnnormalisedReferenceAccumulator;
        private final RatioAccumulator mTumorAccumulator;

        public CobaltAccumulator(final int windowSize, final GenomeRegion region)
        {
            mWindow = new Window(windowSize);
            mRegion = region;

            mReferenceAccumulator = new RatioAccumulator();
            mUnnormalisedReferenceAccumulator = new RatioAccumulator();
            mTumorAccumulator = new RatioAccumulator();
        }

        double referenceMeanRatio()
        {
            return mReferenceAccumulator.meanRatio();
        }

        double unnormalisedReferenceMeanRatio()
        {
            return mUnnormalisedReferenceAccumulator.meanRatio();
        }

        double tumorMeanRatio() { return mTumorAccumulator.meanRatio(); }
        double tumorMedianRatio() { return mTumorAccumulator.medianRatio(); }

        int tumorCount() { return mTumorAccumulator.count(); }

        @Override
        public void accept(final CobaltRatio ratio)
        {
            if(mWindow.end(ratio.position()) <= mRegion.end())
            {
                if(ratio.referenceGCDiploidRatio() < 0)
                    return;

                mReferenceAccumulator.add(ratio.referenceGCDiploidRatio());
                mUnnormalisedReferenceAccumulator.add(ratio.referenceGCRatio());
                mTumorAccumulator.add(ratio.tumorGCRatio(), true);
            }
        }
    }

    static private class RatioAccumulator
    {
        private double mSumRatio;
        private int mCount;
        private final List<Double> mRatios = Lists.newArrayList();

        private double meanRatio()
        {
            return mCount > 0 ? mSumRatio / mCount : 0;
        }

        private double medianRatio()
        {
            if(mRatios.isEmpty())
                return 0;

            int medianIndex = mRatios.size() / 2;

            if((mRatios.size() % 2) == 0)
                return (mRatios.get(medianIndex - 1) + mRatios.get(medianIndex)) * 0.5;
            else
                return mRatios.get(medianIndex);
        }

        private int count()
        {
            return mCount;
        }

        public void add(double ratio)
        {
            add(ratio, false);
        }

        public void add(double ratio, boolean keepValues)
        {
            if(!Doubles.greaterThan(ratio, -1))
                return;

            mCount++;
            mSumRatio += ratio;

            if(keepValues)
            {
                int index = 0;

                while(index < mRatios.size())
                {
                    if(ratio < mRatios.get(index))
                        break;

                    ++index;
                }

                mRatios.add(index, ratio);
            }
        }
    }
}
