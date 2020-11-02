package com.hartwig.hmftools.protect.variants;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariant;
import com.hartwig.hmftools.protect.variants.germline.DriverGermlineVariant;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingEntry;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingModel;
import com.hartwig.hmftools.protect.variants.somatic.DriverSomaticVariant;

import org.jetbrains.annotations.NotNull;

public final class ReportableVariantFactory {

    private ReportableVariantFactory() {
    }

    @NotNull
    public static List<ReportableVariant> reportableGermlineVariants(@NotNull List<ReportableGermlineVariant> variants,
            @NotNull Set<String> genesWithSomaticInactivationEvent, @NotNull GermlineReportingModel germlineReportingModel) {
        List<ReportableVariant> reportableVariants = Lists.newArrayList();

        for (ReportableGermlineVariant variant : variants) {
            GermlineReportingEntry reportingEntry = germlineReportingModel.entryForGene(variant.gene());
            if (reportingEntry != null && isPresentInTumor(variant)) {
                boolean includeVariant;
                String specificHgvsProtein = reportingEntry.reportableHgvsProtein();
                if (specificHgvsProtein != null) {
                    includeVariant = variant.hgvsProtein().equals(specificHgvsProtein);
                } else if (reportingEntry.reportBiallelicOnly()) {
                    boolean filterBiallelic = variant.biallelic();

                    boolean filterGermlineVariantInSameGene = false;
                    for (ReportableGermlineVariant otherVariant : variants) {
                        if (variant != otherVariant && otherVariant.gene().equals(variant.gene())) {
                            filterGermlineVariantInSameGene = true;
                        }
                    }

                    boolean filterSomaticVariantInSameGene = genesWithSomaticInactivationEvent.contains(variant.gene());

                    includeVariant = filterBiallelic || filterGermlineVariantInSameGene || filterSomaticVariantInSameGene;
                } else {
                    includeVariant = true;
                }

                if (includeVariant) {
                    reportableVariants.add(fromGermlineVariant(variant).driverLikelihood(1D).build());
                }
            }
        }

        return reportableVariants;
    }

    @NotNull
    public static List<ReportableVariant> reportableSomaticVariants(@NotNull List<SomaticVariant> variants,
            @NotNull List<DriverCatalog> driverCatalog) {
        Map<String, DriverCatalog> driverCatalogMap = driverCatalog.stream().collect(Collectors.toMap(DriverCatalog::gene, x -> x));

        List<ReportableVariant> result = Lists.newArrayList();
        for (SomaticVariant variant : variants) {
            if (variant.reported()) {
                DriverCatalog geneDriver = driverCatalogMap.get(variant.gene());
                ReportableVariant reportable = fromSomaticVariant(variant).driverLikelihood(geneDriver.driverLikelihood()).build();
                result.add(reportable);
            }
        }

        return result;
    }

    @NotNull
    public static List<ReportableVariant> mergeVariantLists(@NotNull List<ReportableVariant> list1,
            @NotNull List<ReportableVariant> list2) {
        List<ReportableVariant> result = Lists.newArrayList();

        Map<String, Double> maxLikelihood = Maps.newHashMap();
        for (ReportableVariant variant : list1) {
            maxLikelihood.merge(variant.gene(), variant.driverLikelihood(), Math::max);
        }

        for (ReportableVariant variant : list2) {
            maxLikelihood.merge(variant.gene(), variant.driverLikelihood(), Math::max);
        }

        for (ReportableVariant variant : list1) {
            result.add(ImmutableReportableVariant.builder().from(variant).driverLikelihood(maxLikelihood.get(variant.gene())).build());
        }

        for (ReportableVariant variant : list2) {
            result.add(ImmutableReportableVariant.builder().from(variant).driverLikelihood(maxLikelihood.get(variant.gene())).build());
        }

        return result;
    }

    @NotNull
    static List<ReportableVariant> mergeVariantLists(@NotNull List<DriverSomaticVariant> somaticVariantsReport,
            @NotNull List<DriverGermlineVariant> germlineVariantsToReport, @NotNull GermlineReportingModel germlineReportingModel) {
        List<ReportableVariant> allReportableVariants = Lists.newArrayList();
        for (DriverSomaticVariant somaticDriverVariant : somaticVariantsReport) {
            double adjustedDriverLikelihood = somaticDriverVariant.driverLikelihood();
            for (DriverGermlineVariant germlineVariant : germlineVariantsToReport) {
                if (germlineVariant.variant().gene().equals(somaticDriverVariant.variant().gene())) {
                    adjustedDriverLikelihood = Math.max(adjustedDriverLikelihood, germlineVariant.driverLikelihood());
                }
            }

            allReportableVariants.add(fromSomaticVariant(somaticDriverVariant.variant()).driverLikelihood(adjustedDriverLikelihood)
                    .build());
        }

        for (DriverGermlineVariant driverGermlineVariant : germlineVariantsToReport) {
            double adjustedDriverLikelihood = driverGermlineVariant.driverLikelihood();
            for (DriverSomaticVariant somaticVariant : somaticVariantsReport) {
                if (somaticVariant.variant().gene().equals(driverGermlineVariant.variant().gene())) {
                    adjustedDriverLikelihood = Math.max(adjustedDriverLikelihood, somaticVariant.driverLikelihood());
                }
            }

            allReportableVariants.add(fromGermlineVariant(driverGermlineVariant.variant()).driverLikelihood(adjustedDriverLikelihood)
                    .build());
        }

        return allReportableVariants;
    }

    private static boolean isPresentInTumor(@NotNull ReportableGermlineVariant germlineVariant) {
        return calcAlleleCopyNumber(germlineVariant.adjustedCopyNumber(), germlineVariant.adjustedVaf()) >= 0.5;
    }

    @NotNull
    private static ImmutableReportableVariant.Builder fromGermlineVariant(@NotNull ReportableGermlineVariant variant) {
        return ImmutableReportableVariant.builder()
                .source(ReportableVariantSource.GERMLINE)
                .gene(variant.gene())
                .chromosome(variant.chromosome())
                .position(variant.position())
                .ref(variant.ref())
                .alt(variant.alt())
                .canonicalCodingEffect(variant.codingEffect())
                .canonicalHgvsCodingImpact(variant.hgvsCoding())
                .canonicalHgvsProteinImpact(variant.hgvsProtein())
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
                .totalCopyNumber(variant.adjustedCopyNumber())
                .alleleCopyNumber(calcAlleleCopyNumber(variant.adjustedCopyNumber(), variant.adjustedVaf()))
                .hotspot(Hotspot.NON_HOTSPOT)
                .clonalLikelihood(1D)
                .biallelic(variant.biallelic());
    }

    @NotNull
    private static ImmutableReportableVariant.Builder fromSomaticVariant(@NotNull SomaticVariant variant) {
        return ImmutableReportableVariant.builder()
                .source(ReportableVariantSource.SOMATIC)
                .gene(variant.gene())
                .chromosome(variant.chromosome())
                .position(variant.position())
                .ref(variant.ref())
                .alt(variant.alt())
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
