package com.hartwig.hmftools.purple.fitting;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.Doubles.positiveOrZero;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_ANEUPLOIDIC_RATIO_CUTOFF;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_ANEUPLOIDIC_REGION_CUTOFF;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_ANEUPLOIDIC_REGION_MIN_BAF_COUNT;
import static com.hartwig.hmftools.purple.PurpleConstants.HIGHLY_ANEUPLOIDIC_RATIO_CUTOFF;
import static com.hartwig.hmftools.purple.PurpleConstants.HIGHLY_ANEUPLOIDIC_REGION_CUTOFF;
import static com.hartwig.hmftools.purple.PurpleConstants.HIGHLY_ANEUPLOIDIC_REGION_MIN_BAF_COUNT;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.PerChromosomeData;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.purple.region.FittingRegion;

public class AneuploidyDetector
{
    static boolean useRegionToDetectAneuploidy(final PerChromosomeData cobaltChromosomes, final FittingRegion region)
    {
        if(region.bafCount() <= 0)
        {
            return false;
        }

        if(!positiveOrZero(region.observedTumorRatio()))
        {
            return false;
        }

        if(region.germlineStatus() != GermlineStatus.DIPLOID)
        {
            return false;
        }

        return cobaltChromosomes.hasChromosome(region.chromosome());
    }

    private final List<? extends FittingRegion> Regions;
    private final ListMultimap<Chromosome, AmberBAF> AmberData;

    public AneuploidyDetector(
            final List<? extends FittingRegion> regions,
            final ListMultimap<Chromosome, AmberBAF> amberData,
            final PerChromosomeData cobaltChromosomes)
    {
        Regions = regions.stream().filter(region -> useRegionToDetectAneuploidy(cobaltChromosomes, region)).toList();
        AmberData = amberData;
    }

    public boolean hasAneuploidy()
    {
        int totalBafCount = 0;
        int highBafCount = 0;
        for(FittingRegion region : Regions)
        {
            totalBafCount += region.bafCount();
            PPL_LOGGER.debug(format("AR: region %s:%d-%d has %d total points", region.chromosome(), region.start(), region.end(), region.bafCount()));
            if(showsAneuploidy(region))
            {
                PPL_LOGGER.debug(format("AR: region contributes ratio: %d", region.bafCount()));
                highBafCount += region.bafCount();
            }
        }
        PPL_LOGGER.debug(format("AR: high diploid baf count: %d, total diploid baf count: %d", highBafCount, totalBafCount));
        double ratio = (double) highBafCount / totalBafCount;
        PPL_LOGGER.debug(format("aneuploidy ratio: %.3f", ratio));
        if(Doubles.greaterOrEqual(ratio, SOMATIC_FIT_ANEUPLOIDIC_RATIO_CUTOFF))
        {
            return true;
        }
        List<FittingRegion> highlyAneuploidicRegions = new ArrayList<>();
        AmberPointsProvider amberProvider = new AmberPointsProvider(AmberData);
        double requiredPoints = Math.max(HIGHLY_ANEUPLOIDIC_REGION_MIN_BAF_COUNT, (totalBafCount * HIGHLY_ANEUPLOIDIC_RATIO_CUTOFF));
        for(FittingRegion region : Regions)
        {
            if(region.observedBAF() < HIGHLY_ANEUPLOIDIC_REGION_CUTOFF)
            {
                continue;
            }
            if(region.bafCount() < requiredPoints)
            {
                continue;
            }
            List<AmberBAF> amberPoints = amberProvider.nextBlock(region);
            if(isHighlyAneuploidic(amberPoints))
            {
                highlyAneuploidicRegions.add(region);
                PPL_LOGGER.debug(format("Highly aneuploidic region found: %.3f, %d", region.observedBAF(), region.bafCount()));
            }
        }
        return !highlyAneuploidicRegions.isEmpty();
    }

    private boolean isHighlyAneuploidic(final List<AmberBAF> amberPoints)
    {
        Double minVaf = amberPoints.stream().map(AmberBAF::tumorBAF).min(Double::compareTo).orElseThrow();
        return minVaf < 0.3;
    }

    private static boolean showsAneuploidy(final FittingRegion region)
    {
        return region.observedBAF() > SOMATIC_FIT_ANEUPLOIDIC_REGION_CUTOFF
                && region.bafCount() > SOMATIC_FIT_ANEUPLOIDIC_REGION_MIN_BAF_COUNT;
    }
}
