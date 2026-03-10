package com.hartwig.hmftools.purple.fitting;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.Doubles.positiveOrZero;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_ANEUPLOIDIC_RATIO_CUTOFF;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_ANEUPLOIDIC_REGION_CUTOFF;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_ANEUPLOIDIC_REGION_MIN_BAF_COUNT;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.util.List;

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

    private final List<FittingRegion> Regions;

    public AneuploidyDetector(final List<FittingRegion> regions, final PerChromosomeData cobaltChromosomes)
    {
        Regions = regions.stream().filter(region -> useRegionToDetectAneuploidy(cobaltChromosomes, region)).toList();
    }

    public boolean hasAneuploidy()
    {
        int totalBafCount = 0;
        int highBafCount = 0;
        for(FittingRegion region : Regions)
        {
            totalBafCount += region.bafCount();
            PPL_LOGGER.debug(format("AR: region %s:%d-%d has %d total points", region.chromosome(), region.start(), region.end(), region.bafCount()));
            if(region.observedBAF() > SOMATIC_FIT_ANEUPLOIDIC_REGION_CUTOFF
                    && region.bafCount() > SOMATIC_FIT_ANEUPLOIDIC_REGION_MIN_BAF_COUNT)
            {
                PPL_LOGGER.debug(format("AR: region contributes ratio: %d", region.bafCount()));
                highBafCount += region.bafCount();
            }
        }
        PPL_LOGGER.debug(format("AR: high diploid baf count: %d, total diploid baf count: %d", highBafCount, totalBafCount));
        double ratio = (double) highBafCount / totalBafCount;
        PPL_LOGGER.debug(format("aneuploidy ratio: %3f", ratio));
        return Doubles.greaterOrEqual(ratio, SOMATIC_FIT_ANEUPLOIDIC_RATIO_CUTOFF);
    }
}
