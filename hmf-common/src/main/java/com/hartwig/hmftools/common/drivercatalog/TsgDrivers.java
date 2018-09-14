package com.hartwig.hmftools.common.drivercatalog;

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
    public static List<DriverCatalog> tsgDrivers(@NotNull final Map<String, DndsDriverImpactLikelihood> likelihoodsByGene,
            @NotNull final List<EnrichedSomaticVariant> variants) {

        final List<DriverCatalog> driverCatalog = Lists.newArrayList();

        final Map<String, List<EnrichedSomaticVariant>> codingVariants =
                DriverCatalogFactory.codingVariantsByGene(likelihoodsByGene.keySet(), variants);

        return driverCatalog;
    }

    @NotNull
    static DriverCatalog geneDriver(long sampleSNVCount, @NotNull final DndsDriverGeneLikelihood likelihood,
            @NotNull final List<EnrichedSomaticVariant> codingVariants) {
        codingVariants.sort(new TsgImpactComparator());

        final ImmutableDriverCatalog.Builder builder = ImmutableDriverCatalog.builder()
                .gene(likelihood.gene())
                .category(DriverCategory.TSG)
                .driverLikelihood(1)
                .dndsLikelihood(0)
                .driver(codingVariants.size() > 1 ? DriverType.MULTI_HIT : DriverType.SINGLE_HIT);

        if (codingVariants.stream().anyMatch(SomaticVariant::hotspot)) {
            return builder.driver(DriverType.HOTSPOT).build();
        }

        if (codingVariants.stream().anyMatch(PurityAdjustedSomaticVariant::biallelic)) {
            return builder.driver(DriverType.BIALLELIC).build();
        }

        final DndsDriverImpactLikelihood firstImpactLikelihood = impactLikelihood(likelihood, codingVariants.get(0));
        if (codingVariants.size() == 1) {
            // SingleHit
            return builder.dndsLikelihood(firstImpactLikelihood.dndsLikelihood())
                    .driverLikelihood(singleHit(sampleSNVCount, firstImpactLikelihood))
                    .build();
        }

        // Must be multi hit
        final DndsDriverImpactLikelihood secondImpactLikelihood = impactLikelihood(likelihood, codingVariants.get(1));

        return builder.dndsLikelihood(Math.max(firstImpactLikelihood.dndsLikelihood(), secondImpactLikelihood.dndsLikelihood()))
                .driverLikelihood(DriverCatalogFactory.probabilityDriverVariant(sampleSNVCount,
                        firstImpactLikelihood,
                        secondImpactLikelihood))
                .build();
    }

    static DndsDriverImpactLikelihood impactLikelihood(@NotNull final DndsDriverGeneLikelihood dndsLikelihood,
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

    static double singleHit(long sampleSNVCount, @NotNull final DndsDriverImpactLikelihood likelihood) {

        return Doubles.positive(likelihood.dndsLikelihood())
                ? DriverCatalogFactory.probabilityDriverVariant(sampleSNVCount, likelihood)
                : 0;
    }
}
