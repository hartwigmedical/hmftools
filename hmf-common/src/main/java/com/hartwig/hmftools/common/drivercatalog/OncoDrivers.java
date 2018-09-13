package com.hartwig.hmftools.common.drivercatalog;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.dnds.DndsDriverLikelihood;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.jetbrains.annotations.NotNull;

public class OncoDrivers {

    static final int MAX_REPEAT_COUNT = 7;

    @NotNull
    public static List<DriverCatalog> oncoDrivers(@NotNull final Map<String, DndsDriverLikelihood> likelihoodsByGene,
            @NotNull final List<EnrichedSomaticVariant> variants) {

        final List<DriverCatalog> driverCatalog = Lists.newArrayList();

        final Map<VariantType, Long> variantTypeCounts =
                variants.stream().filter(x -> !x.isFiltered()).collect(Collectors.groupingBy(SomaticVariant::type, Collectors.counting()));
        long sampleSNVCount = variantTypeCounts.getOrDefault(VariantType.SNP, 0L);

        final Map<String, List<EnrichedSomaticVariant>> variantsByGene = variants.stream()
                .filter(x -> likelihoodsByGene.keySet().contains(x.gene()))
                .filter(x -> x.canonicalCodingEffect() != CodingEffect.SYNONYMOUS && x.canonicalCodingEffect() != CodingEffect.NONE)
                .collect(Collectors.groupingBy(SomaticVariant::gene));

        for (String gene : variantsByGene.keySet()) {

            final DndsDriverLikelihood geneLikelihood = likelihoodsByGene.get(gene);
            final List<EnrichedSomaticVariant> geneVariants = variantsByGene.get(gene);

            driverCatalog.add(geneDriver(sampleSNVCount, geneLikelihood, geneVariants));
        }

        return driverCatalog;
    }

    @NotNull
    static DriverCatalog geneDriver(long sampleSNVCount, @NotNull final DndsDriverLikelihood likelihood,
            @NotNull final List<EnrichedSomaticVariant> geneVariants) {

        long missenseVariants = geneVariants.stream().filter(OncoDrivers::isMissense).count();

        final ImmutableDriverCatalog.Builder builder = ImmutableDriverCatalog.builder()
                .gene(likelihood.gene())
                .category(DriverCategory.ONCO)
                .driverLikelihood(1)
                .dndsLikelihood(missenseVariants > 0 ? likelihood.missenseUnadjustedDriverLikelihood() : 0)
                .driver(missenseVariants > 1 ? DriverType.MULTI_HIT : DriverType.SINGLE_HIT);

        if (geneVariants.stream().anyMatch(SomaticVariant::hotspot)) {
            return builder.driver(DriverType.HOTSPOT).build();
        }

        // TODO: Add NEAR_HOTSPOT HERE

        if (geneVariants.stream().anyMatch(OncoDrivers::isInframeIndel)) {
            return builder.driver(DriverType.INFRAME).build();
        }

        final double driverLikelihood = Doubles.positive(likelihood.missenseUnadjustedDriverLikelihood()) && missenseVariants > 0
                ? missenseProbabilityDriverVariant(sampleSNVCount, likelihood)
                : 0;

        return builder.driverLikelihood(driverLikelihood).build();
    }

    private static boolean isInframeIndel(@NotNull final EnrichedSomaticVariant variant) {
        return variant.type() == VariantType.INDEL && variant.canonicalCodingEffect() == CodingEffect.MISSENSE
                && variant.repeatCount() <= MAX_REPEAT_COUNT;
    }

    private static boolean isMissense(@NotNull final EnrichedSomaticVariant variant) {
        return variant.type() != VariantType.INDEL && variant.canonicalCodingEffect() == CodingEffect.MISSENSE;
    }

    public static double missenseProbabilityDriverVariant(long sampleSNVCount, @NotNull final DndsDriverLikelihood likelihood) {

        double lambda = sampleSNVCount * likelihood.missenseProbabilityVariantNonDriverFactor();
        PoissonDistribution poissonDistribution = new PoissonDistribution(lambda);

        double pDriver = likelihood.missenseProbabilityDriver();
        double pVariantNonDriver = 1 - poissonDistribution.cumulativeProbability(0);

        return pDriver / (pDriver + pVariantNonDriver * (1 - pDriver));
    }

}
