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

public class DriverCatalogFactory {

    public static <T extends SomaticVariant> Map<DriverImpact, Long> driverImpactCount(@NotNull final List<T> variants) {
        return variants.stream().collect(Collectors.groupingBy(DriverImpact::select, Collectors.counting()));
    }

    public static <T extends SomaticVariant> Map<VariantType, Long> variantTypeCount(@NotNull final List<T> variants) {
        return variants.stream().collect(Collectors.groupingBy(SomaticVariant::type, Collectors.counting()));
    }

    static double probabilityDriverVariant(long sampleSNVCount, @NotNull final DndsDriverImpactLikelihood likelihood) {
        double lambda = sampleSNVCount * likelihood.pVariantNonDriverFactor();
        if (Doubles.isZero(lambda)) {
            return likelihood.dndsLikelihood();
        }

        PoissonDistribution poissonDistribution = new PoissonDistribution(lambda);

        double pVariantNonDriver = 1 - poissonDistribution.cumulativeProbability(0);
        return likelihood.pDriver() / (likelihood.pDriver() + pVariantNonDriver * (1 - likelihood.pDriver()));
    }

    static double probabilityDriverVariant(long firstVariantTypeCount, long secondVariantTypeCount,
            @NotNull final DndsDriverImpactLikelihood firstLikelihood, @NotNull final DndsDriverImpactLikelihood secondLikelihood) {

        final double pDriver = Math.max(firstLikelihood.pDriver(), secondLikelihood.pDriver());

        final double pVariantNonDriver;
        if (firstLikelihood.equals(secondLikelihood)) {
            double lambda = firstVariantTypeCount * firstLikelihood.pVariantNonDriverFactor();
            PoissonDistribution poissonDistribution = new PoissonDistribution(lambda);
            pVariantNonDriver = 1 - poissonDistribution.cumulativeProbability(1);
        } else {
            double pVariantNonDriver1 =
                    1 - new PoissonDistribution(firstVariantTypeCount * firstLikelihood.pVariantNonDriverFactor()).cumulativeProbability(0);
            double pVariantNonDriver2 = 1 - new PoissonDistribution(
                    secondVariantTypeCount * secondLikelihood.pVariantNonDriverFactor()).cumulativeProbability(0);
            pVariantNonDriver = pVariantNonDriver1 * pVariantNonDriver2;
        }

        return pDriver / (pDriver + pVariantNonDriver * (1 - pDriver));
    }

}
