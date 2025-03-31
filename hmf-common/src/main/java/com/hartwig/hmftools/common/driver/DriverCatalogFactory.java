package com.hartwig.hmftools.common.driver;

import com.hartwig.hmftools.common.driver.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.math3.distribution.PoissonDistribution;

public final class DriverCatalogFactory
{
    private DriverCatalogFactory() {}

    public static double probabilityDriverVariant(long sampleSNVCount, final DndsDriverImpactLikelihood likelihood)
    {
        return probabilityDriverVariantSameImpact(0, sampleSNVCount, likelihood);
    }

    private static double probabilityDriverVariantSameImpact(int count, long sampleSNVCount, final DndsDriverImpactLikelihood likelihood)
    {
        double lambda = sampleSNVCount * likelihood.passengersPerMutation();
        if(Doubles.isZero(lambda))
            return 0.0;

        PoissonDistribution poissonDistribution = new PoissonDistribution(lambda);

        double pVariantNonDriver = 1 - poissonDistribution.cumulativeProbability(count);
        return likelihood.driversPerSample() / (likelihood.driversPerSample() + pVariantNonDriver * (1 - likelihood.driversPerSample()));
    }

    public static double probabilityDriverVariant(
            long firstVariantTypeCount, long secondVariantTypeCount,
            final DndsDriverImpactLikelihood firstLikelihood, final DndsDriverImpactLikelihood secondLikelihood)
    {
        if(firstLikelihood.equals(secondLikelihood))
        {
            return probabilityDriverVariantSameImpact(1, firstVariantTypeCount, firstLikelihood);
        }

        double lambda1 = firstVariantTypeCount * firstLikelihood.passengersPerMutation();
        double lambda2 = secondVariantTypeCount * secondLikelihood.passengersPerMutation();
        if(Doubles.isZero(lambda1) || Doubles.isZero(lambda2))
        {
            return Math.max(probabilityDriverVariant(firstVariantTypeCount, firstLikelihood),
                    probabilityDriverVariant(secondVariantTypeCount, secondLikelihood));
        }

        final double pDriver = Math.max(firstLikelihood.driversPerSample(), secondLikelihood.driversPerSample());
        final double pVariantNonDriver1 = 1 - new PoissonDistribution(lambda1).cumulativeProbability(0);
        final double pVariantNonDriver2 = 1 - new PoissonDistribution(lambda2).cumulativeProbability(0);
        final double pVariantNonDriver = pVariantNonDriver1 * pVariantNonDriver2;

        return pDriver / (pDriver + pVariantNonDriver * (1 - pDriver));
    }

    public static DriverCatalog createCopyNumberDriver(
            DriverCategory category, DriverType driver, final LikelihoodMethod likelihoodMethod, final boolean biallelic,
            final GeneCopyNumber geneCopyNumber)
    {
        return ImmutableDriverCatalog.builder()
                .chromosome(geneCopyNumber.chromosome())
                .chromosomeBand(geneCopyNumber.chromosomeBand())
                .gene(geneCopyNumber.geneName())
                .transcript(geneCopyNumber.transName())
                .isCanonical(geneCopyNumber.isCanonical())
                .missense(0)
                .nonsense(0)
                .inframe(0)
                .frameshift(0)
                .splice(0)
                .driverLikelihood(1)
                .driver(driver)
                .likelihoodMethod(likelihoodMethod)
                .category(category)
                .biallelic(biallelic)
                .minCopyNumber(geneCopyNumber.minCopyNumber())
                .maxCopyNumber(geneCopyNumber.maxCopyNumber())
                .build();
    }

}
