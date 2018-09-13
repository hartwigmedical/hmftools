package com.hartwig.hmftools.common.drivercatalog;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.jetbrains.annotations.NotNull;

public class DriverCatalogFactory {

    public static <T extends SomaticVariant> Map<String, List<T>> codingVariantsByGene(@NotNull final Set<String> genes,
            @NotNull final List<T> variants) {

        return variants.stream()
                .filter(x -> genes.contains(x.gene()))
                .filter(x -> x.canonicalCodingEffect() != CodingEffect.SYNONYMOUS && x.canonicalCodingEffect() != CodingEffect.NONE)
                .collect(Collectors.groupingBy(SomaticVariant::gene));
    }

    static double probabilityDriverVariant(boolean multiHit, long sampleSNVCount, @NotNull final DndsDriverImpactLikelihood likelihood) {
        double lambda = sampleSNVCount * likelihood.pVariantNonDriverFactor();
        PoissonDistribution poissonDistribution = new PoissonDistribution(lambda);

        double pVariantNonDriver = 1 - poissonDistribution.cumulativeProbability(multiHit ? 1 : 0);
        return likelihood.pDriver() / (likelihood.pDriver() + pVariantNonDriver * (1 - likelihood.pDriver()));
    }

}
