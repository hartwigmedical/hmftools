package com.hartwig.hmftools.common.drivercatalog;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.genome.region.TranscriptRegion;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class OncoDrivers {

    static final int MAX_REPEAT_COUNT = 7;

    private OncoDrivers() {
    }

    @NotNull
    static List<DriverCatalog> drivers(@NotNull final Map<String, DndsDriverImpactLikelihood> likelihoodsByGene,
            @NotNull final List<SomaticVariant> variants, @NotNull final List<GeneCopyNumber> geneCopyNumberList,
            final Map<VariantType, Long> variantTypeCounts) {
        final List<DriverCatalog> driverCatalog = Lists.newArrayList();

        final Map<String, GeneCopyNumber> geneCopyNumbers =
                geneCopyNumberList.stream().collect(Collectors.toMap(TranscriptRegion::gene, x -> x));
        long sampleSNVCount = variantTypeCounts.getOrDefault(VariantType.SNP, 0L);

        final Map<String, List<SomaticVariant>> codingVariants = oncogenicVariantsByGene(likelihoodsByGene.keySet(), variants);

        for (String gene : codingVariants.keySet()) {
            final DndsDriverImpactLikelihood geneMissenseLikelihood = likelihoodsByGene.get(gene);
            final List<SomaticVariant> geneVariants = codingVariants.get(gene);

            driverCatalog.add(geneDriver(sampleSNVCount, gene, geneMissenseLikelihood, geneVariants, geneCopyNumbers.get(gene)));
        }

        return driverCatalog;
    }

    @NotNull
    private static Map<String, List<SomaticVariant>> oncogenicVariantsByGene(@NotNull final Set<String> genes,
            @NotNull final List<SomaticVariant> variants) {
        final Predicate<SomaticVariant> oncoPredicate = oncoVariant(genes);
        return variants.stream().filter(oncoPredicate).collect(Collectors.groupingBy(SomaticVariant::gene));
    }

    @NotNull
    static Predicate<SomaticVariant> oncoVariant(@NotNull final Set<String> genes) {
        return x -> genes.contains(x.gene()) && (isMissense(x) || isInframeIndel(x) || x.hotspot() == Hotspot.HOTSPOT);
    }

    @NotNull
    static DriverCatalog geneDriver(long sampleSNVCount, @NotNull final String gene,
            @NotNull final DndsDriverImpactLikelihood missenseLikelihood, @NotNull final List<SomaticVariant> geneVariants,
            @Nullable GeneCopyNumber geneCopyNumber) {
        final Map<DriverImpact, Long> variantCounts = DriverCatalogFactory.driverImpactCount(geneVariants);
        long missenseVariants = variantCounts.getOrDefault(DriverImpact.MISSENSE, 0L);
        long nonsenseVariants = variantCounts.getOrDefault(DriverImpact.NONSENSE, 0L);
        long spliceVariants = variantCounts.getOrDefault(DriverImpact.SPLICE, 0L);
        long inframeVariants = variantCounts.getOrDefault(DriverImpact.INFRAME, 0L);
        long frameshiftVariants = variantCounts.getOrDefault(DriverImpact.FRAMESHIFT, 0L);

        final ImmutableDriverCatalog.Builder builder = ImmutableDriverCatalog.builder()
                .chromosome(geneVariants.get(0).chromosome())
                .chromosomeBand(geneCopyNumber == null ? Strings.EMPTY : geneCopyNumber.chromosomeBand())
                .gene(gene)
                .driver(DriverType.MUTATION)
                .category(DriverCategory.ONCO)
                .driverLikelihood(1)
                .dndsLikelihood(missenseVariants > 0 ? missenseLikelihood.dndsLikelihood() : 0)
                .missense(missenseVariants)
                .nonsense(nonsenseVariants)
                .splice(spliceVariants)
                .inframe(inframeVariants)
                .frameshift(frameshiftVariants)
                .biallelic(geneVariants.stream().anyMatch(SomaticVariant::biallelic))
                .minCopyNumber(geneCopyNumber == null ? 0 : geneCopyNumber.minCopyNumber())
                .maxCopyNumber(geneCopyNumber == null ? 0 : geneCopyNumber.maxCopyNumber())
                .likelihoodMethod(missenseVariants > 0 ? LikelihoodMethod.DNDS : LikelihoodMethod.NONE);

        if (geneVariants.stream().anyMatch(SomaticVariant::isHotspot)) {
            return builder.likelihoodMethod(LikelihoodMethod.HOTSPOT).build();
        }

        if (geneVariants.stream().anyMatch(OncoDrivers::isInframeIndel)) {
            return builder.likelihoodMethod(LikelihoodMethod.INFRAME).build();
        }

        final double driverLikelihood =
                Doubles.positive(missenseLikelihood.dndsLikelihood()) && missenseVariants > 0 ? missenseProbabilityDriverVariant(
                        sampleSNVCount,
                        missenseLikelihood) : 0;

        return builder.driverLikelihood(driverLikelihood).build();
    }

    private static boolean isMissense(@NotNull final SomaticVariant variant) {
        return variant.type() != VariantType.INDEL && variant.canonicalCodingEffect() == CodingEffect.MISSENSE;
    }

    private static boolean isInframeIndel(@NotNull final SomaticVariant variant) {
        return variant.type() == VariantType.INDEL && variant.canonicalCodingEffect() == CodingEffect.MISSENSE
                && variant.repeatCount() <= MAX_REPEAT_COUNT;
    }

    private static double missenseProbabilityDriverVariant(long sampleSNVCount, @NotNull final DndsDriverImpactLikelihood likelihood) {
        return DriverCatalogFactory.probabilityDriverVariant(sampleSNVCount, likelihood);
    }
}
