package com.hartwig.hmftools.common.variant;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogKey;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogMap;
import com.hartwig.hmftools.common.drivercatalog.DriverType;

import org.jetbrains.annotations.NotNull;

public final class ReportableVariantFactory {

    private ReportableVariantFactory() {
    }

    @NotNull
    public static List<ReportableVariant> toReportableGermlineVariants(@NotNull List<SomaticVariant> variants,
            @NotNull List<DriverCatalog> germlineDriverCatalog) {
        List<DriverCatalog> germlineMutationCatalog =
                germlineDriverCatalog.stream().filter(x -> x.driver() == DriverType.GERMLINE_MUTATION).collect(Collectors.toList());
        return toReportableVariants(variants, germlineMutationCatalog, ReportableVariantSource.GERMLINE);
    }

    @NotNull
    public static List<ReportableVariant> toReportableSomaticVariants(@NotNull List<SomaticVariant> variants,
            @NotNull List<DriverCatalog> somaticDriverCatalog) {
        List<DriverCatalog> somaticMutationCatalog =
                somaticDriverCatalog.stream().filter(x -> x.driver() == DriverType.MUTATION).collect(Collectors.toList());
        return toReportableVariants(variants, somaticMutationCatalog, ReportableVariantSource.SOMATIC);
    }

    @NotNull
    private static List<ReportableVariant> toReportableVariants(@NotNull List<SomaticVariant> variants,
            @NotNull List<DriverCatalog> driverCatalog, @NotNull ReportableVariantSource source) {
        Map<DriverCatalogKey, DriverCatalog> geneDriverMap = DriverCatalogMap.toDriverMap(driverCatalog);
        List<ReportableVariant> result = Lists.newArrayList();

        Set<DriverCatalogKey> keys = DriverCatalogKey.buildUniqueKeysSet(geneDriverMap);

        for (SomaticVariant variant : variants) {
            if (variant.reported()) {
                ImmutableReportableVariant.Builder reportable = ImmutableReportableVariant.builder();
                for (DriverCatalogKey key : keys) {
                    DriverCatalog geneDriver = geneDriverMap.get(key);
                    if (geneDriver == null) {
                        throw new IllegalStateException("Could not find driver entry for variant on gene '" + variant.gene() + "'");
                    }
                    SomaticVariant variantCorrect = ImmutableSomaticVariantImpl.builder().from(variant).gene(geneDriver.gene()).build();
                    ImmutableReportableVariant.Builder build =
                            fromVariant(variantCorrect, source).driverLikelihood(geneDriver.driverLikelihood())
                                    .transcript(geneDriver.transcript())
                                    .isCanonical(geneDriver.isCanonical());
                    reportable.from(build.build());
                }

                result.add(reportable.build());
            }
        }

        return result;
    }

    @NotNull
    public static List<ReportableVariant> mergeVariantLists(@NotNull List<ReportableVariant> list1,
            @NotNull List<ReportableVariant> list2) {
        List<ReportableVariant> result = Lists.newArrayList();

        Map<String, Double> maxLikelihoodPerGene = Maps.newHashMap();
        for (ReportableVariant variant : list1) {
            maxLikelihoodPerGene.merge(variant.gene(), variant.driverLikelihood(), Math::max);
        }

        for (ReportableVariant variant : list2) {
            maxLikelihoodPerGene.merge(variant.gene(), variant.driverLikelihood(), Math::max);
        }

        for (ReportableVariant variant : list1) {
            result.add(ImmutableReportableVariant.builder()
                    .from(variant)
                    .driverLikelihood(maxLikelihoodPerGene.get(variant.gene()))
                    .build());
        }

        for (ReportableVariant variant : list2) {
            result.add(ImmutableReportableVariant.builder()
                    .from(variant)
                    .driverLikelihood(maxLikelihoodPerGene.get(variant.gene()))
                    .build());
        }

        return result;
    }

    @NotNull
    public static ImmutableReportableVariant.Builder fromVariant(@NotNull SomaticVariant variant, @NotNull ReportableVariantSource source) {
        return ImmutableReportableVariant.builder()
                .type(variant.type())
                .source(source)
                .gene(variant.gene())
                .chromosome(variant.chromosome())
                .position(variant.position())
                .ref(variant.ref())
                .alt(variant.alt())
                .otherReportedEffects(variant.otherReportedEffects())
                .canonicalTranscript(variant.canonicalTranscript())
                .canonicalEffect(variant.canonicalEffect())
                .canonicalCodingEffect(variant.canonicalCodingEffect())
                .canonicalHgvsCodingImpact(variant.canonicalHgvsCodingImpact())
                .canonicalHgvsProteinImpact(variant.canonicalHgvsProteinImpact())
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
                .totalCopyNumber(variant.adjustedCopyNumber())
                .minorAlleleCopyNumber(variant.minorAlleleCopyNumber())
                .alleleCopyNumber(calcAlleleCopyNumber(variant.adjustedCopyNumber(), variant.adjustedVAF()))
                .hotspot(variant.hotspot())
                .clonalLikelihood(variant.clonalLikelihood())
                .biallelic(variant.biallelic())
                .genotypeStatus(variant.genotypeStatus())
                .localPhaseSet(variant.topLocalPhaseSet());
    }

    private static double calcAlleleCopyNumber(double adjustedCopyNumber, double adjustedVAF) {
        return adjustedCopyNumber * Math.max(0, Math.min(1, adjustedVAF));
    }
}
