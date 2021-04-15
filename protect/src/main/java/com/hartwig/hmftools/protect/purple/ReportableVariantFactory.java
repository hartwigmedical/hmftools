package com.hartwig.hmftools.protect.purple;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.protect.germline.GermlineReportingEntry;
import com.hartwig.hmftools.protect.germline.GermlineReportingModel;

import org.jetbrains.annotations.NotNull;

public final class ReportableVariantFactory {

    private ReportableVariantFactory() {
    }

    @NotNull
    public static List<ReportableVariant> reportableGermlineVariants(@NotNull List<SomaticVariant> variants,
            @NotNull List<DriverCatalog> germlineDriverCatalog, @NotNull GermlineReportingModel germlineReportingModel) {
        return reportableVariants(variants, germlineDriverCatalog, ReportableVariantSource.GERMLINE, germlineReportingModel);
    }

    @NotNull
    public static List<ReportableVariant> reportableSomaticVariants(@NotNull List<SomaticVariant> variants,
            @NotNull List<DriverCatalog> somaticDriverCatalog, @NotNull GermlineReportingModel germlineReportingModel) {
        List<DriverCatalog> somaticMutationCatalog =
                somaticDriverCatalog.stream().filter(x -> x.driver() == DriverType.MUTATION).collect(Collectors.toList());
        return reportableVariants(variants, somaticMutationCatalog, ReportableVariantSource.SOMATIC, germlineReportingModel);
    }

    @NotNull
    private static List<ReportableVariant> reportableVariants(@NotNull List<SomaticVariant> variants,
            @NotNull List<DriverCatalog> driverCatalog, @NotNull ReportableVariantSource source, @NotNull GermlineReportingModel germlineReportingModel) {
        Map<String, DriverCatalog> geneDriverMap = driverCatalog.stream().collect(Collectors.toMap(DriverCatalog::gene, x -> x));

        List<ReportableVariant> result = Lists.newArrayList();
        for (SomaticVariant variant : variants) {
            if (variant.reported()) {
                boolean includeVariant = true;
                if (source == ReportableVariantSource.GERMLINE) {
                    GermlineReportingEntry reportingEntry = germlineReportingModel.entryForGene(variant.gene());
                    if (reportingEntry != null) {
                        String exclusiveHgvsProteinFilter = reportingEntry.exclusiveHgvsProteinFilter();
                        if (exclusiveHgvsProteinFilter != null) {
                            includeVariant = variant.canonicalHgvsProteinImpact().equals(exclusiveHgvsProteinFilter);
                        }
                    }
                }

                if (includeVariant) {
                    DriverCatalog geneDriver = geneDriverMap.get(variant.gene());

                    if (geneDriver == null) {
                        throw new IllegalStateException("Could not find driver entry for variant on gene '" + variant.gene() + "'");
                    }

                    ReportableVariant reportable = fromVariant(variant, source).driverLikelihood(geneDriver.driverLikelihood()).build();
                    result.add(reportable);
                }
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
    private static ImmutableReportableVariant.Builder fromVariant(@NotNull SomaticVariant variant,
            @NotNull ReportableVariantSource source) {
        return ImmutableReportableVariant.builder()
                .type(variant.type())
                .source(source)
                .gene(variant.gene())
                .chromosome(variant.chromosome())
                .position(variant.position())
                .ref(variant.ref())
                .alt(variant.alt())
                .canonicalTranscript(variant.canonicalTranscript())
                .canonicalCodingEffect(variant.canonicalCodingEffect())
                .canonicalHgvsCodingImpact(variant.canonicalHgvsCodingImpact())
                .canonicalHgvsProteinImpact(variant.canonicalHgvsProteinImpact())
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
                .totalCopyNumber(variant.adjustedCopyNumber())
                .alleleCopyNumber(calcAlleleCopyNumber(variant.adjustedCopyNumber(), variant.adjustedVAF()))
                .hotspot(variant.hotspot())
                .clonalLikelihood(variant.clonalLikelihood())
                .biallelic(variant.biallelic());
    }

    private static double calcAlleleCopyNumber(double adjustedCopyNumber, double adjustedVAF) {
        return adjustedCopyNumber * Math.max(0, Math.min(1, adjustedVAF));
    }
}
