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

public class TsgDrivers {

    private static final int MAX_REPEAT_COUNT = 7;

    @NotNull
    public static List<DriverCatalog> tsgDrivers(@NotNull final Map<String, DndsDriverLikelihood> likelihoodsByGene,
            @NotNull final List<EnrichedSomaticVariant> variants) {

        final List<DriverCatalog> driverCatalog = Lists.newArrayList();

        final Map<VariantType, Long> variantTypeCounts =
                variants.stream().filter(x -> !x.isFiltered()).collect(Collectors.groupingBy(SomaticVariant::type, Collectors.counting()));
        long sampleSNVCount = variantTypeCounts.getOrDefault(VariantType.SNP, 0L);

        final Map<String, List<EnrichedSomaticVariant>> missenseVariantsByGene = variants.stream()
                .filter(x -> likelihoodsByGene.keySet().contains(x.gene()))
                .filter(x -> x.canonicalCodingEffect().equals(CodingEffect.MISSENSE))
                .collect(Collectors.groupingBy(SomaticVariant::gene));

        for (String gene : missenseVariantsByGene.keySet()) {
            final DndsDriverLikelihood likelihood = likelihoodsByGene.get(gene);
            final ImmutableDriverCatalog.Builder builder = ImmutableDriverCatalog.builder()
                    .gene(gene)
                    .category(DriverCategory.TSG)
                    .driverLikelihood(1)
                    .dndsLikelihood(likelihood.missenseUnadjustedDriverLikelihood());

            final List<EnrichedSomaticVariant> geneVariants = missenseVariantsByGene.get(gene);

            boolean isHotspot = geneVariants.stream().anyMatch(SomaticVariant::hotspot);
            boolean isBiallelic = geneVariants.stream().anyMatch(EnrichedSomaticVariant::biallelic);
            boolean isNonHotspotNorBiallelic = geneVariants.stream().anyMatch(x -> !x.hotspot() && !x.biallelic());

            if (isHotspot) {
                driverCatalog.add(builder.driver(DriverType.HOTSPOT).build());
            }

            if (isBiallelic) {
                driverCatalog.add(builder.driver(DriverType.BIALLELIC).build());
            }

            if (isNonHotspotNorBiallelic) {

                if (geneVariants.size() == 1) {
                    // SingleHit

                } else {
                    // MultiHit so sort and look at worst two.

                }

            }

            /// RUBBISH UNDER HERE

            final List<EnrichedSomaticVariant> inframe = geneVariants.stream()
                    .filter(x -> x.type() == VariantType.INDEL && x.repeatCount() <= MAX_REPEAT_COUNT)
                    .collect(Collectors.toList());
            if (!inframe.isEmpty()) {
                driverCatalog.add(builder.driver(DriverType.INFRAME).build());
            }

            final List<EnrichedSomaticVariant> pointMutations =
                    geneVariants.stream().filter(x -> x.type() != VariantType.INDEL).collect(Collectors.toList());

            if (Doubles.greaterThan(likelihood.missenseUnadjustedDriverLikelihood(), 0) && pointMutations.stream()
                    .anyMatch(x -> !x.hotspot())) {
                final double missenseDriverLikelihood = missenseProbabilityDriverVariant(sampleSNVCount, likelihood);
                driverCatalog.add(builder.driver(geneVariants.size() > 1 ? DriverType.MULTI_HIT : DriverType.SINGLE_HIT)
                        .driverLikelihood(missenseDriverLikelihood)
                        .build());
            }
        }

        return driverCatalog;
    }

    public static double missenseProbabilityDriverVariant(long sampleSNVCount, @NotNull final DndsDriverLikelihood likelihood) {

        double lambda = sampleSNVCount * likelihood.missenseProbabilityVariantNonDriverFactor();
        PoissonDistribution poissonDistribution = new PoissonDistribution(lambda);

        double pDriver = likelihood.missenseProbabilityDriver();
        double pVariantNonDriver = 1 - poissonDistribution.cumulativeProbability(0);

        return pDriver / (pDriver + pVariantNonDriver * (1 - pDriver));
    }

}
