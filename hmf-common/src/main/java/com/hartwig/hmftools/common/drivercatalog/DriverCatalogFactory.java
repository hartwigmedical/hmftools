package com.hartwig.hmftools.common.drivercatalog;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.jetbrains.annotations.NotNull;

public final class DriverCatalogFactory {

    private DriverCatalogFactory() {
    }

    @NotNull
    public static <T extends SomaticVariant> Map<DriverImpact, Long> driverImpactCount(@NotNull final List<T> variants) {
        return variants.stream().collect(Collectors.groupingBy(DriverImpact::select, Collectors.counting()));
    }

    @NotNull
    static <T extends SomaticVariant> Map<VariantType, Long> variantTypeCount(@NotNull final List<T> variants) {
        return variants.stream().collect(Collectors.groupingBy(SomaticVariant::type, Collectors.counting()));
    }

    static double probabilityDriverVariant(long sampleSNVCount, @NotNull final DndsDriverImpactLikelihood likelihood) {
        return probabilityDriverVariantSameImpact(0, sampleSNVCount, likelihood);
    }

    private static double probabilityDriverVariantSameImpact(int count, long sampleSNVCount,
            @NotNull final DndsDriverImpactLikelihood likelihood) {
        double lambda = sampleSNVCount * likelihood.pVariantNonDriverFactor();
        if (Doubles.isZero(lambda)) {
            return likelihood.dndsLikelihood();
        }

        PoissonDistribution poissonDistribution = new PoissonDistribution(lambda);

        double pVariantNonDriver = 1 - poissonDistribution.cumulativeProbability(count);
        return likelihood.pDriver() / (likelihood.pDriver() + pVariantNonDriver * (1 - likelihood.pDriver()));
    }

    static double probabilityDriverVariant(long firstVariantTypeCount, long secondVariantTypeCount,
            @NotNull final DndsDriverImpactLikelihood firstLikelihood, @NotNull final DndsDriverImpactLikelihood secondLikelihood) {
        if (firstLikelihood.equals(secondLikelihood)) {
            return probabilityDriverVariantSameImpact(1, firstVariantTypeCount, firstLikelihood);
        }

        double lambda1 = firstVariantTypeCount * firstLikelihood.pVariantNonDriverFactor();
        double lambda2 = secondVariantTypeCount * secondLikelihood.pVariantNonDriverFactor();
        if (Doubles.isZero(lambda1) || Doubles.isZero(lambda2)) {
            return Math.max(firstLikelihood.dndsLikelihood(), secondLikelihood.dndsLikelihood());
        }

        final double pDriver = Math.max(firstLikelihood.pDriver(), secondLikelihood.pDriver());
        final double pVariantNonDriver1 = 1 - new PoissonDistribution(lambda1).cumulativeProbability(0);
        final double pVariantNonDriver2 = 1 - new PoissonDistribution(lambda2).cumulativeProbability(0);
        final double pVariantNonDriver = pVariantNonDriver1 * pVariantNonDriver2;

        return pDriver / (pDriver + pVariantNonDriver * (1 - pDriver));
    }
}
