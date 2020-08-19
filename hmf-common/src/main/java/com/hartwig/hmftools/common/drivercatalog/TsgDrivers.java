package com.hartwig.hmftools.common.drivercatalog;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.genome.region.TranscriptRegion;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class TsgDrivers {

    private final ReportablePredicate reportablePredicate;
    private final Map<String, DndsDriverGeneLikelihood> likelihoodsByGene;

    public TsgDrivers(DriverGenePanel genePanel) {
        likelihoodsByGene = genePanel.tsgLikelihood();
        reportablePredicate = new ReportablePredicate(DriverCategory.TSG, genePanel);
    }

    @NotNull
    List<DriverCatalog> drivers(@NotNull final List<SomaticVariant> variants, @NotNull final List<GeneCopyNumber> geneCopyNumberList,
            final Map<VariantType, Long> variantTypeCounts, final Map<VariantType, Long> variantTypeCountsBiallelic,
            final Map<VariantType, Long> variantTypeCountsNonBiallelic) {
        final Map<String, GeneCopyNumber> geneCopyNumbers =
                geneCopyNumberList.stream().collect(Collectors.toMap(TranscriptRegion::gene, x -> x));
        final List<DriverCatalog> driverCatalog = Lists.newArrayList();

        final Map<String, List<SomaticVariant>> codingVariants = codingVariantsByGene(likelihoodsByGene.keySet(), variants);

        for (String gene : codingVariants.keySet()) {
            final DndsDriverGeneLikelihood likelihood = likelihoodsByGene.get(gene);

            final List<SomaticVariant> geneVariants = codingVariants.get(gene);
            driverCatalog.add(geneDriver(likelihood,
                    geneVariants,
                    variantTypeCounts,
                    variantTypeCountsBiallelic,
                    variantTypeCountsNonBiallelic,
                    geneCopyNumbers.get(gene)));
        }

        return driverCatalog;
    }

    @NotNull
    static DriverCatalog geneDriver(@NotNull final DndsDriverGeneLikelihood likelihood, @NotNull final List<SomaticVariant> geneVariants,
            @NotNull final Map<VariantType, Long> standardCounts, @NotNull final Map<VariantType, Long> biallelicCounts,
            @NotNull final Map<VariantType, Long> nonBiallelicCounts, @Nullable GeneCopyNumber geneCopyNumber) {
        geneVariants.sort(new TsgImpactComparator());

        final Map<DriverImpact, Long> variantCounts = DriverCatalogFactory.driverImpactCount(geneVariants);
        long missenseVariants = variantCounts.getOrDefault(DriverImpact.MISSENSE, 0L);
        long nonsenseVariants = variantCounts.getOrDefault(DriverImpact.NONSENSE, 0L);
        long spliceVariants = variantCounts.getOrDefault(DriverImpact.SPLICE, 0L);
        long inframeVariants = variantCounts.getOrDefault(DriverImpact.INFRAME, 0L);
        long frameshiftVariants = variantCounts.getOrDefault(DriverImpact.FRAMESHIFT, 0L);

        final ImmutableDriverCatalog.Builder builder = ImmutableDriverCatalog.builder()
                .chromosome(geneVariants.get(0).chromosome())
                .chromosomeBand(geneCopyNumber == null ? Strings.EMPTY : geneCopyNumber.chromosomeBand())
                .gene(likelihood.gene())
                .driver(DriverType.MUTATION)
                .category(DriverCategory.TSG)
                .driverLikelihood(1)
                .dndsLikelihood(0.0)
                .missense(missenseVariants)
                .nonsense(nonsenseVariants)
                .splice(spliceVariants)
                .inframe(inframeVariants)
                .frameshift(frameshiftVariants)
                .biallelic(geneVariants.stream().anyMatch(SomaticVariant::biallelic))
                .minCopyNumber(geneCopyNumber == null ? 0 : geneCopyNumber.minCopyNumber())
                .maxCopyNumber(geneCopyNumber == null ? 0 : geneCopyNumber.maxCopyNumber())
                .likelihoodMethod(LikelihoodMethod.DNDS);

        if (geneVariants.stream().anyMatch(SomaticVariant::isHotspot)) {
            return builder.likelihoodMethod(LikelihoodMethod.HOTSPOT).build();
        }

        if (geneVariants.stream().anyMatch(x -> x.biallelic() && !DriverImpact.isMissense(x))) {
            return builder.likelihoodMethod(LikelihoodMethod.BIALLELIC).build();
        }

        final DriverImpact firstImpact = DriverImpact.select(geneVariants.get(0));
        final DndsDriverImpactLikelihood firstImpactLikelihood = likelihood.select(geneVariants.get(0).biallelic(), firstImpact);
        final long firstVariantTypeCount =
                variantCount(likelihood.useBiallelic(), geneVariants.get(0), standardCounts, biallelicCounts, nonBiallelicCounts);

        if (geneVariants.size() == 1) {
            return builder.driverLikelihood(singleHit(firstVariantTypeCount, firstImpactLikelihood)).build();
        }

        // MultiHit
        final DriverImpact secondImpact = DriverImpact.select(geneVariants.get(1));
        final DndsDriverImpactLikelihood secondImpactLikelihood = likelihood.select(geneVariants.get(1).biallelic(), secondImpact);
        final long secondVariantTypeCount =
                variantCount(likelihood.useBiallelic(), geneVariants.get(1), standardCounts, biallelicCounts, nonBiallelicCounts);

        return builder.driverLikelihood(DriverCatalogFactory.probabilityDriverVariant(firstVariantTypeCount,
                secondVariantTypeCount,
                firstImpactLikelihood,
                secondImpactLikelihood)).build();
    }

    private static double singleHit(long sampleCount, @NotNull final DndsDriverImpactLikelihood likelihood) {
        return DriverCatalogFactory.probabilityDriverVariant(sampleCount, likelihood);
    }

    @NotNull
    private <T extends SomaticVariant> Map<String, List<T>> codingVariantsByGene(@NotNull final Set<String> genes,
            @NotNull final List<T> variants) {
        return variants.stream().filter(reportablePredicate).collect(Collectors.groupingBy(SomaticVariant::gene));
    }

    private static long variantCount(boolean useBiallelic, @NotNull final SomaticVariant variant,
            @NotNull final Map<VariantType, Long> standard, @NotNull final Map<VariantType, Long> biallelic,
            @NotNull final Map<VariantType, Long> nonBiallelic) {
        final Map<VariantType, Long> map;
        if (!useBiallelic) {
            map = standard;
        } else if (variant.biallelic()) {
            map = biallelic;
        } else {
            map = nonBiallelic;
        }

        return variant.type() == VariantType.INDEL ? map.getOrDefault(VariantType.INDEL, 0L) : map.getOrDefault(VariantType.SNP, 0L);
    }

}
