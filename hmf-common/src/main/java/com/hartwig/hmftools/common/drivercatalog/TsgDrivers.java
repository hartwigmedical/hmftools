package com.hartwig.hmftools.common.drivercatalog;

import static com.hartwig.hmftools.common.drivercatalog.DriverCatalogFactory.variantTypeCount;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;

public class TsgDrivers {

    @NotNull
    public static List<DriverCatalog> tsgDrivers(@NotNull final Map<String, DndsDriverGeneLikelihood> likelihoodsByGene,
            @NotNull final List<EnrichedSomaticVariant> variants) {

        final List<DriverCatalog> driverCatalog = Lists.newArrayList();

        final Map<VariantType, Long> variantTypeCounts = variantTypeCount(variants);
        long sampleSNVCount = variantTypeCounts.getOrDefault(VariantType.SNP, 0L);
        long sampleIndelCount = variantTypeCounts.getOrDefault(VariantType.INDEL, 0L);

        final Map<String, List<EnrichedSomaticVariant>> codingVariants =
                DriverCatalogFactory.codingVariantsByGene(likelihoodsByGene.keySet(), variants);

        for (String gene : codingVariants.keySet()) {
            final List<EnrichedSomaticVariant> geneVariants = codingVariants.get(gene);
            driverCatalog.add(geneDriver(sampleSNVCount, sampleIndelCount, likelihoodsByGene.get(gene), geneVariants));
        }

        return driverCatalog;
    }

    @NotNull
    static DriverCatalog geneDriver(long sampleSNVCount, long sampleIndelCount, @NotNull final DndsDriverGeneLikelihood likelihood,
            @NotNull final List<EnrichedSomaticVariant> codingVariants) {
        codingVariants.sort(new TsgImpactComparator());

        final Map<DriverImpact, Long> variantCounts = DriverCatalogFactory.driverImpactCount(codingVariants);
        long missenseVariants = variantCounts.getOrDefault(DriverImpact.MISSENSE, 0L);
        long nonsenseVariants = variantCounts.getOrDefault(DriverImpact.NONSENSE, 0L);
        long spliceVariants = variantCounts.getOrDefault(DriverImpact.SPLICE, 0L);
        long inframeVariants = variantCounts.getOrDefault(DriverImpact.INFRAME, 0L);
        long frameshiftVariants = variantCounts.getOrDefault(DriverImpact.FRAMESHIFT, 0L);

        final double maxDndsLikelihood = codingVariants.stream()
                .map(x -> impactLikelihood(likelihood, x))
                .mapToDouble(DndsDriverImpactLikelihood::dndsLikelihood)
                .max()
                .orElse(0);

        final ImmutableDriverCatalog.Builder builder = ImmutableDriverCatalog.builder()
                .gene(likelihood.gene())
                .category(DriverCategory.TSG)
                .driverLikelihood(1)
                .dndsLikelihood(maxDndsLikelihood)
                .missense(missenseVariants)
                .nonsense(nonsenseVariants)
                .splice(spliceVariants)
                .inframe(inframeVariants)
                .frameshift(frameshiftVariants)
                .driver(DriverType.DNDS);

        if (codingVariants.stream().anyMatch(SomaticVariant::hotspot)) {
            return builder.driver(DriverType.HOTSPOT).build();
        }

        if (codingVariants.stream().anyMatch(PurityAdjustedSomaticVariant::biallelic)) {
            return builder.driver(DriverType.BIALLELIC).build();
        }

        final DndsDriverImpactLikelihood firstImpactLikelihood = impactLikelihood(likelihood, codingVariants.get(0));
        final long firstVariantTypeCount = codingVariants.get(0).type() == VariantType.INDEL ? sampleIndelCount : sampleSNVCount;

        if (codingVariants.size() == 1) {
            return builder.dndsLikelihood(firstImpactLikelihood.dndsLikelihood())
                    .driverLikelihood(singleHit(firstVariantTypeCount, firstImpactLikelihood))
                    .build();
        }

        // MultiHit
        final DndsDriverImpactLikelihood secondImpactLikelihood = impactLikelihood(likelihood, codingVariants.get(1));
        final long secondVariantTypeCount = codingVariants.get(1).type() == VariantType.INDEL ? sampleIndelCount : sampleSNVCount;

        return builder.dndsLikelihood(Math.max(firstImpactLikelihood.dndsLikelihood(), secondImpactLikelihood.dndsLikelihood()))
                .driverLikelihood(DriverCatalogFactory.probabilityDriverVariant(firstVariantTypeCount,
                        secondVariantTypeCount,
                        firstImpactLikelihood,
                        secondImpactLikelihood))
                .build();
    }

    private static DndsDriverImpactLikelihood impactLikelihood(@NotNull final DndsDriverGeneLikelihood dndsLikelihood,
            @NotNull final SomaticVariant variant) {

        if (variant.type() == VariantType.INDEL) {
            return dndsLikelihood.indel();
        } else if (variant.canonicalCodingEffect() == CodingEffect.MISSENSE) {
            return dndsLikelihood.missense();
        } else if (variant.canonicalCodingEffect() == CodingEffect.NONSENSE_OR_FRAMESHIFT) {
            return dndsLikelihood.nonsense();
        }
        return dndsLikelihood.splice();
    }

    private static double singleHit(long sampleCount, @NotNull final DndsDriverImpactLikelihood likelihood) {

        return Doubles.positive(likelihood.dndsLikelihood()) ? DriverCatalogFactory.probabilityDriverVariant(sampleCount, likelihood) : 0;
    }
}
