package com.hartwig.hmftools.purple.region;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.purple.GermlineStatus.CENTROMETIC;
import static com.hartwig.hmftools.common.purple.GermlineStatus.DIPLOID;
import static com.hartwig.hmftools.common.purple.GermlineStatus.EXCLUDED;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleConstants.CENTROMERIC_WIDTH;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_AMP_DEL_EXCLUSION_CHR_1;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_AMP_DEL_EXCLUSION_CHR_17;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_AMP_DEL_EXCLUSION_CHR_19;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_AMP_DEL_EXCLUSION_CHR_9;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_AMP_RATIO;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_DEL_MIN_LENGTH;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_DEL_RATIO;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.genome.region.Window;
import com.hartwig.hmftools.common.immune.ImmuneRegions;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.purple.segment.PurpleSegment;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;

public class ObservedRegionFactory
{
    private final int mWindowSize;
    private final CobaltChromosomes mCobaltChromosomes;
    private final GermlineStatusCalcs mStatusFactory;

    private static List<ChrBaseRegion> EXCLUDED_IMMUNE_REGIONS = Lists.newArrayList();
    private static List<ChrBaseRegion> CENTROMETRIC_REGIONS = Lists.newArrayList();

    private static List<ChrBaseRegion> GERMLINE_AMP_DEL_EXCLUSIONS = Lists.newArrayList();

    public ObservedRegionFactory(final int windowSize, final CobaltChromosomes cobaltChromosomes)
    {
        mWindowSize = windowSize;
        mCobaltChromosomes = cobaltChromosomes;
        mStatusFactory = new GermlineStatusCalcs(cobaltChromosomes);
    }

    public static void setSpecificRegions(
            final RefGenomeVersion refGenomeVersion, final Map<Chromosome,GenomePosition> chrLengths,
            final Map<Chromosome,GenomePosition> centromeres)
    {
        EXCLUDED_IMMUNE_REGIONS.addAll(ImmuneRegions.getIgRegions(refGenomeVersion));
        EXCLUDED_IMMUNE_REGIONS.addAll(ImmuneRegions.getTrRegions(refGenomeVersion));

        int halfWidth = CENTROMERIC_WIDTH / 2;

        for(GenomePosition centromere : centromeres.values())
        {
            int centromereStart = centromere.position() - halfWidth;
            int centromereEnd = centromere.position() + halfWidth;
            CENTROMETRIC_REGIONS.add(new ChrBaseRegion(centromere.chromosome(), centromereStart, centromereEnd));
        }

        GERMLINE_AMP_DEL_EXCLUSIONS.add(new ChrBaseRegion(
                refGenomeVersion.versionedChromosome("1"), 1, GERMLINE_AMP_DEL_EXCLUSION_CHR_1));

        GERMLINE_AMP_DEL_EXCLUSIONS.add(new ChrBaseRegion(
                refGenomeVersion.versionedChromosome("9"), GERMLINE_AMP_DEL_EXCLUSION_CHR_9,
                chrLengths.get(HumanChromosome._9).position()));

        GERMLINE_AMP_DEL_EXCLUSIONS.add(new ChrBaseRegion(
                refGenomeVersion.versionedChromosome("17"), GERMLINE_AMP_DEL_EXCLUSION_CHR_17,
                chrLengths.get(HumanChromosome._17).position()));

        GERMLINE_AMP_DEL_EXCLUSIONS.add(new ChrBaseRegion(
                refGenomeVersion.versionedChromosome("19"), 1, GERMLINE_AMP_DEL_EXCLUSION_CHR_19));
    }

    public List<ObservedRegion> formObservedRegions(
            final List<PurpleSegment> regions, final Multimap<Chromosome, AmberBAF> bafs,
            final Map<Chromosome,List<CobaltRatio>> ratios, final Multimap<Chromosome,GCProfile> gcProfiles)
    {
        List<ObservedRegion> observedRegions = Lists.newArrayList();

        GenomePositionSelector<CobaltRatio> cobaltSelector = GenomePositionSelectorFactory.create(ratios);
        GenomePositionSelector<AmberBAF> bafSelector = GenomePositionSelectorFactory.create(bafs);
        GenomeRegionSelector<GCProfile> gcSelector = GenomeRegionSelectorFactory.createImproved(gcProfiles);

        List<Integer> candidateGermlineAmpDelRegions = Lists.newArrayList();

        for(PurpleSegment region : regions)
        {
            final BAFAccumulator baf = new BAFAccumulator();
            final CobaltAccumulator cobalt = new CobaltAccumulator(mWindowSize, region);
            final GCAccumulator gc = new GCAccumulator(region);

            bafSelector.select(region, baf);
            cobaltSelector.select(region, cobalt);
            gcSelector.select(region, gc);

            double tumorRatio = cobalt.tumorMedianRatio();
            double normalRatio = cobalt.referenceMeanRatio();
            int depthWindowCount = cobalt.tumorCount();

            GermlineStatus germlineStatus = getGermlineStatus(region, normalRatio, tumorRatio, depthWindowCount);

            ObservedRegion observedRegion = new ObservedRegion(
                    region.chromosome(), region.start(), region.end(), region.RatioSupport, region.Support, baf.count(), baf.medianBaf(),
                    depthWindowCount, tumorRatio, normalRatio, cobalt.unnormalisedReferenceMeanRatio(), germlineStatus,
                    region.SvCluster, gc.averageGCContent(), region.MinStart, region.MaxStart);

            if(observedRegion.start() > observedRegion.end()
            || !positionsWithin(region.MinStart, region.MaxStart, observedRegion.start(), observedRegion.end()))
            {
                PPL_LOGGER.error("invalid observed region: {}", observedRegion);
            }

            observedRegions.add(observedRegion);

            double rawNormalRatio = cobalt.unnormalisedReferenceMeanRatio();

            if(isGermlineAmpDelCandidate(region, germlineStatus, rawNormalRatio, normalRatio))
            {
                // record index of the candidate region
                candidateGermlineAmpDelRegions.add(observedRegions.size() - 1);
            }
        }

        findGermlineAmpDelRegions(observedRegions, candidateGermlineAmpDelRegions);

        extendMinSupport(observedRegions);
        return observedRegions;
    }

    private GermlineStatus getGermlineStatus(final PurpleSegment region, double normalRatio, double tumorRatio, int depthWindowCount)
    {
        if(EXCLUDED_IMMUNE_REGIONS.stream()
                .anyMatch(x -> x.Chromosome.equals(region.Chromosome) && positionsWithin(region.start(), region.end(), x.start(), x.end())))
        {
            return GermlineStatus.EXCLUDED;
        }

        if(CENTROMETRIC_REGIONS.stream()
                .anyMatch(x -> x.Chromosome.equals(region.Chromosome) && positionsWithin(region.start(), region.end(), x.start(), x.end())))
        {
            return GermlineStatus.CENTROMETIC;
        }

        return mStatusFactory.calcStatus(region.chromosome(), normalRatio, tumorRatio, depthWindowCount);
    }

    private boolean isGermlineAmpDelCandidate(
            final PurpleSegment region, final GermlineStatus germlineStatus, double rawNormalRatio, double normalRatio)
    {
        if(germlineStatus != DIPLOID || normalRatio <= 0 || HumanChromosome.fromString(region.Chromosome).isAllosome())
            return false;

        if(GERMLINE_AMP_DEL_EXCLUSIONS.stream()
                .anyMatch(x -> x.Chromosome.equals(region.Chromosome) && positionsWithin(region.start(), region.end(), x.start(), x.end())))
        {
            return false;
        }

        return isGermlineAmpDelRatio(rawNormalRatio, normalRatio);
    }

    private static boolean isGermlineAmpDelCandidate(final ObservedRegion region)
    {
        return isGermlineAmpDelRatio(region.unnormalisedObservedNormalRatio(), region.observedNormalRatio());
    }

    private static boolean isGermlineAmpDelRatio(double rawNormalRatio, double normalRatio)
    {
        if(normalRatio == 0)
            return false;

        double ratio = rawNormalRatio / normalRatio;
        return ratio < GERMLINE_DEL_RATIO || ratio > GERMLINE_AMP_RATIO;
    }

    private void findGermlineAmpDelRegions(final List<ObservedRegion> observedRegions, final List<Integer> candidateGermlineAmpDelRegions)
    {
        // test for diploid regions with evidence of germline AMPs or DELs
        if(candidateGermlineAmpDelRegions.isEmpty())
            return;

        Set<Integer> processedCandidateIndices = Sets.newHashSet();

        for(Integer candidateRegionIndex : candidateGermlineAmpDelRegions)
        {
            if(processedCandidateIndices.contains(candidateRegionIndex))
                continue;

            ObservedRegion candidateRegion = observedRegions.get(candidateRegionIndex);

            int candidateRangeMin = 0;
            int candidateRangeMax = 0;
            int candidateIndexMin = 0;
            int candidateIndexMax = 0;

            // walk forwards and backwards to the nearest diploid regions or end of the arm
            for(int j = candidateRegionIndex - 1; j >= 0; --j)
            {
                ObservedRegion nextRegion = observedRegions.get(j);

                if(!nextRegion.chromosome().equals(candidateRegion.chromosome()))
                    break;

                if(nextRegion.support() == SegmentSupport.CENTROMERE
                || nextRegion.germlineStatus() == DIPLOID && !isGermlineAmpDelCandidate(nextRegion))
                {
                    candidateRangeMin = nextRegion.end();
                    break;
                }
                else if(nextRegion.support() == SegmentSupport.TELOMERE)
                {
                    candidateRangeMin = nextRegion.start();
                    candidateIndexMin = j;
                    break;
                }

                candidateIndexMin = j;
            }

            for(int j = candidateRegionIndex + 1; j < observedRegions.size(); ++j)
            {
                ObservedRegion nextRegion = observedRegions.get(j);

                if(!nextRegion.chromosome().equals(candidateRegion.chromosome()))
                {
                    // telomere segment status hasn't been set, so take the previous
                    candidateIndexMax = j - 1;
                    candidateRangeMax = observedRegions.get(candidateIndexMax).end();

                    break;
                }

                if(nextRegion.support() == SegmentSupport.CENTROMERE
                || nextRegion.germlineStatus() == DIPLOID && !isGermlineAmpDelCandidate(nextRegion))
                {
                    candidateRangeMax = nextRegion.start();
                    break;
                }
                else if(nextRegion.support() == SegmentSupport.TELOMERE)
                {
                    candidateRangeMax = nextRegion.end();
                    candidateIndexMax = j;
                    break;
                }

                candidateIndexMax = j;
            }

            if(candidateRangeMax - candidateRangeMin >= GERMLINE_DEL_MIN_LENGTH)
            {
                PPL_LOGGER.info("germline event from region({}:{}-{}) normalRatio({}) range({}-{})",
                        candidateRegion.chromosome(), candidateRegion.start(), candidateRegion.end(),
                        format("%.2f unnorm=%.2f", candidateRegion.observedNormalRatio(),
                                candidateRegion.unnormalisedObservedNormalRatio()),
                        candidateRangeMin, candidateRangeMax);

                // override the normalised ratio for all regions within these bounds
                for(int j = candidateIndexMin; j <= candidateIndexMax; ++j)
                {
                    ObservedRegion region = observedRegions.get(j);
                    region.setObservedNormalRatio(region.unnormalisedObservedNormalRatio());

                    // recalc germline status
                    if(region.germlineStatus() != EXCLUDED && region.germlineStatus() != CENTROMETIC)
                    {
                        GermlineStatus newGermlineStatus = mStatusFactory.calcStatus(
                                region.chromosome(), region.observedNormalRatio(), region.observedTumorRatio(), region.depthWindowCount());

                        region.setGermlineStatus(newGermlineStatus);
                    }

                    processedCandidateIndices.add(j);
                }
            }
        }
    }

    public static void extendMinSupport(final List<ObservedRegion> observedRegions)
    {
        for(int i = 0; i < observedRegions.size(); i++)
        {
            final ObservedRegion target = observedRegions.get(i);
            if(target.support() == SegmentSupport.NONE && target.germlineStatus() == DIPLOID)
            {
                for(int j = i - 1; j >= 0; j--)
                {
                    final ObservedRegion prior = observedRegions.get(j);
                    if(prior.germlineStatus() == DIPLOID)
                    {
                        break;
                    }

                    target.setMinStart(min(target.minStart(), prior.start()));

                    if(prior.support() != SegmentSupport.NONE)
                    {
                        break;
                    }
                }
            }
        }
    }

    private class BAFAccumulator implements Consumer<AmberBAF>
    {
        private int mCount;
        final private List<Double> mBafs = Lists.newArrayList();

        @Override
        public void accept(final AmberBAF baf)
        {
            if(mCobaltChromosomes.hasChromosome(baf.Chromosome))
            {
                CobaltChromosome cobaltChromosome = mCobaltChromosomes.get(baf.Chromosome);
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
