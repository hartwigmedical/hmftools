package com.hartwig.hmftools.protect.variants;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariant;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingEntry;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingModel;

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
                String exclusiveHgvsProteinFilter = reportingEntry.exclusiveHgvsProteinFilter();
                if (exclusiveHgvsProteinFilter != null) {
                    includeVariant = variant.hgvsProtein().equals(exclusiveHgvsProteinFilter);
                } else if (reportingEntry.reportBiallelicOnly()) {
                    boolean filterBiallelic = variant.biallelic();

                    boolean filterGermlineVariantInSameGene = false;
                    for (ReportableGermlineVariant otherVariant : variants) {
                        if (variant != otherVariant && otherVariant.gene().equals(variant.gene())) {
                            filterGermlineVariantInSameGene = true;
                        }
                    }

                    boolean filterSomaticEventInSameGene = genesWithSomaticInactivationEvent.contains(variant.gene());

                    includeVariant = filterBiallelic || filterGermlineVariantInSameGene || filterSomaticEventInSameGene;
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
        Map<String, DriverCatalog> mutationDriverMap = Maps.newHashMap();

        for (DriverCatalog entry : driverCatalog) {
            if (entry.driver() == DriverType.MUTATION) {
                mutationDriverMap.put(entry.gene(), entry);
            }
        }

        List<ReportableVariant> result = Lists.newArrayList();
        for (SomaticVariant variant : variants) {
            if (variant.reported()) {
                DriverCatalog geneDriver = mutationDriverMap.get(variant.gene());
                assert geneDriver != null;

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

    private static boolean isPresentInTumor(@NotNull ReportableGermlineVariant germlineVariant) {
        return calcAlleleCopyNumber(germlineVariant.adjustedCopyNumber(), germlineVariant.adjustedVaf()) >= 0.5;
    }

    private static double calcAlleleCopyNumber(double adjustedCopyNumber, double adjustedVAF) {
        return adjustedCopyNumber * Math.max(0, Math.min(1, adjustedVAF));
    }
}
