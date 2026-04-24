package com.hartwig.hmftools.purple.fitting;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.Doubles.positiveOrZero;
import static com.hartwig.hmftools.purple.PurpleConstants.HIGHLY_ANEUPLOIDIC_REFIT_BAF_CUTOFF;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_ANEUPLOIDIC_RATIO_CUTOFF;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_ANEUPLOIDIC_REGION_CUTOFF;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_ANEUPLOIDIC_REGION_MIN_BAF_COUNT;
import static com.hartwig.hmftools.purple.PurpleConstants.HIGHLY_ANEUPLOIDIC_RATIO_CUTOFF;
import static com.hartwig.hmftools.purple.PurpleConstants.HIGHLY_ANEUPLOIDIC_REGION_CUTOFF;
import static com.hartwig.hmftools.purple.PurpleConstants.HIGHLY_ANEUPLOIDIC_REGION_MIN_BAF_COUNT;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.PerChromosomeData;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.purple.region.FittingRegion;

public class AneuploidyDetector
{
    private final List<? extends FittingRegion> mRegions;
    private final ListMultimap<Chromosome,AmberBAF> mAmberData;
    private boolean mHasHighBafRegion;
    private int mTotalBafCount;

    public AneuploidyDetector(
            final List<? extends FittingRegion> regions, final ListMultimap<Chromosome, AmberBAF> amberData,
            final PerChromosomeData cobaltChromosomes)
    {
        mRegions = regions.stream().filter(region -> useRegionToDetectAneuploidy(cobaltChromosomes, region)).collect(Collectors.toList());
        mAmberData = amberData;
        mHasHighBafRegion = false;
        mTotalBafCount = 0;
    }

    public boolean hasAneuploidy()
    {
        int highBafCount = 0;

        for(FittingRegion region : mRegions)
        {
            mTotalBafCount += region.bafCount();
            // PPL_LOGGER.trace(format("region %s:%d-%d has %d total points", region.chromosome(), region.start(), region.end(), region.bafCount()));

            if(showsAneuploidy(region))
            {
                // PPL_LOGGER.trace(format("region contributes ratio: %d", region.bafCount()));
                highBafCount += region.bafCount();

                if(region.observedBAF() > HIGHLY_ANEUPLOIDIC_REFIT_BAF_CUTOFF)
                    mHasHighBafRegion = true;
            }
        }

        // PPL_LOGGER.trace(format("high diploid baf count: %d, total diploid baf count: %d", highBafCount, mTotalBafCount));

        double ratio = (double)highBafCount / mTotalBafCount;

        PPL_LOGGER.debug(format("aneuploidy ratio(%.3f)", ratio));

        if(Doubles.greaterOrEqual(ratio, SOMATIC_FIT_ANEUPLOIDIC_RATIO_CUTOFF))
            return true;

        List<FittingRegion> highlyAneuploidicRegions = new ArrayList<>();
        AmberPointsProvider amberProvider = new AmberPointsProvider(mAmberData);
        double requiredPoints = Math.max(HIGHLY_ANEUPLOIDIC_REGION_MIN_BAF_COUNT, (mTotalBafCount * HIGHLY_ANEUPLOIDIC_RATIO_CUTOFF));

        for(FittingRegion region : mRegions)
        {
            if(region.observedBAF() < HIGHLY_ANEUPLOIDIC_REGION_CUTOFF)
                continue;

            if(region.bafCount() < requiredPoints)
                continue;

            List<AmberBAF> amberPoints = amberProvider.nextBlock(region);

            if(isHighlyAneuploidic(amberPoints))
            {
                highlyAneuploidicRegions.add(region);

                PPL_LOGGER.debug(format("region(%s) with highly aneuploidic found: baf(observed=%.3f count=%d)",
                        region.chrBaseRegion(), region.observedBAF(), region.bafCount()));
            }
        }

        return !highlyAneuploidicRegions.isEmpty();
    }

    public boolean hasHighBafRegion() { return mHasHighBafRegion; }

    public Double calculateBafPurity()
    {
        if(mRegions.isEmpty())
            return null;

        Collections.sort(mRegions, Comparator.comparingDouble(x -> x.observedBAF()));

        int cumulativeTotal = 0;
        double maxPercThreshold = mTotalBafCount * (1 - HIGHLY_ANEUPLOIDIC_RATIO_CUTOFF);

        for(FittingRegion fittingRegion : mRegions)
        {
            cumulativeTotal += fittingRegion.bafCount();

            if(cumulativeTotal > maxPercThreshold && fittingRegion.bafCount() >= HIGHLY_ANEUPLOIDIC_REGION_MIN_BAF_COUNT)
            {
                if(fittingRegion.observedBAF() < SOMATIC_FIT_ANEUPLOIDIC_REGION_CUTOFF)
                    return null;

                return fittingRegion.observedBAF() * 2 - 1;
            }
        }

        return null;
    }

    private boolean isHighlyAneuploidic(final List<AmberBAF> amberPoints)
    {
        if(amberPoints.isEmpty())
            return false;

        // ensure a point above and below 0.3/0.7
        boolean hasHigh = false;
        boolean hasLow = false;

        for(AmberBAF amberBAF : amberPoints)
        {
            if(amberBAF.TumorBAF > HIGHLY_ANEUPLOIDIC_REGION_CUTOFF)
                hasHigh = true;
            else if(amberBAF.TumorBAF < (1 - HIGHLY_ANEUPLOIDIC_REGION_CUTOFF))
                hasLow = true;

            if(hasHigh && hasLow)
                return true;
        }

        return false;
    }

    private static boolean showsAneuploidy(final FittingRegion region)
    {
        return region.observedBAF() > SOMATIC_FIT_ANEUPLOIDIC_REGION_CUTOFF
            && region.bafCount() > SOMATIC_FIT_ANEUPLOIDIC_REGION_MIN_BAF_COUNT;
    }

    @VisibleForTesting
    public static boolean useRegionToDetectAneuploidy(final PerChromosomeData cobaltChromosomes, final FittingRegion region)
    {
        if(region.bafCount() <= 0)
            return false;

        if(!positiveOrZero(region.observedTumorRatio()))
            return false;

        if(region.germlineStatus() != GermlineStatus.DIPLOID)
            return false;

        return cobaltChromosomes.hasChromosome(region.chromosome());
    }
}
