package com.hartwig.hmftools.purple.fittingsnv;

import static java.lang.Math.floor;
import static java.lang.String.format;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.config.PurpleConstants.SNV_READJUST_NTH_PERC;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class SomaticReadjustmentFit
{

    public SomaticReadjustmentFit()
    {

    }

    public void assessReadjustment(final List<SomaticVariant> variants, final FittedPurity standardSomaticPurity)
    {
        double standardPurity = standardSomaticPurity.purity();

        List<CandidateOutlier> candidateOutliers = Lists.newArrayList();

        for(SomaticVariant variant : variants)
        {
            if(!HumanChromosome.fromString(variant.chromosome()).isAutosome())
                continue;

            double expectedAlleleCount = variant.totalReadCount() * standardPurity * 0.5;

            if(variant.alleleReadCount() > expectedAlleleCount)
            {
                PoissonDistribution poissonDist = new PoissonDistribution(expectedAlleleCount);
                double probability = 1 - poissonDist.cumulativeProbability(variant.alleleReadCount() - 1);

                candidateOutliers.add(new CandidateOutlier(variant, probability));
            }

            PPL_LOGGER.trace(format("hotspot(%s:%d) vaf(%.3f %d/%d)",
                    variant.chromosome(), variant.position(),
                    variant.alleleFrequency(), variant.alleleReadCount(), variant.totalReadCount()));
        }

        Collections.sort(candidateOutliers, Comparator.comparingDouble(x -> x.Probability));

        // find the nth percentile
        int nthPercIndex = (int)floor(candidateOutliers.size() * SNV_READJUST_NTH_PERC);

        if(nthPercIndex >= 1)
        {
            CandidateOutlier nthCandidate = candidateOutliers.get(nthPercIndex);

            PPL_LOGGER.debug("nth outlier candidate({})", nthCandidate.Variant);
        }
    }

    private class CandidateOutlier
    {
        public final SomaticVariant Variant;
        public final double Probability;

        public CandidateOutlier(final SomaticVariant variant, final double probability)
        {
            Variant = variant;
            Probability = probability;
        }

        public String toString() { return format("%s prob(%.6f)", Variant, Probability); }
    }

}
