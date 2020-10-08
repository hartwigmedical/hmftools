package com.hartwig.hmftools.protect.variants.germline;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariant;
import com.hartwig.hmftools.protect.variants.somatic.DriverSomaticVariant;

import org.jetbrains.annotations.NotNull;

public final class FilterGermlineVariants {

    private FilterGermlineVariants() {
    }

    @NotNull
    public static List<DriverGermlineVariant> filterGermlineVariantsForReporting(
            @NotNull List<ReportableGermlineVariant> germlineVariants, @NotNull GermlineReportingModel germlineReportingModel,
            @NotNull List<GeneCopyNumber> allGeneCopyNumbers, @NotNull List<DriverSomaticVariant> somaticDriverVariantsToReport,
            @NotNull ChordStatus chordStatus) {
        List<DriverGermlineVariant> reportableGermlineVariants = Lists.newArrayList();

        Set<String> reportableGermlineGenes = germlineReportingModel.reportableGermlineGenes();
        for (ReportableGermlineVariant germlineVariant : presentInTumorOnly(germlineVariants)) {
            if (reportableGermlineGenes.contains(germlineVariant.gene())) {
                // TODO: This will be cleaned up in upcoming germline variant reporting logic
                if (germlineVariant.gene().equals("KIT")) {
                    reportableGermlineVariants.add(reportableGermlineVariantWithDriverLikelihood(germlineVariant, 1.0));
                } else {
                    // Only report germline variants on TSGs if there is a 2nd hit or CHORD suggests HRD
                    boolean filterBiallelic = germlineVariant.biallelic();

                    boolean filterMinCopyNumberTumor = false;
                    GeneCopyNumber geneCopyNumber = lookupGeneCopyNumber(allGeneCopyNumbers, germlineVariant.gene());
                    if (Math.round(geneCopyNumber.minCopyNumber()) <= 1 && (Math.round(germlineVariant.adjustedCopyNumber()) >= 2)) {
                        filterMinCopyNumberTumor = true;
                    }

                    boolean filterSomaticVariantInSameGene = false;
                    for (DriverSomaticVariant driverVariant : somaticDriverVariantsToReport) {
                        if (driverVariant.variant().gene().equals(germlineVariant.gene())) {
                            filterSomaticVariantInSameGene = true;
                        }
                    }

                    boolean filterGermlineVariantInSameGene = false;
                    for (ReportableGermlineVariant variant : germlineVariants) {
                        if (variant != germlineVariant && variant.gene().equals(germlineVariant.gene())) {
                            filterGermlineVariantInSameGene = true;
                        }
                    }

                    if (filterBiallelic || filterSomaticVariantInSameGene || filterGermlineVariantInSameGene) {
                        reportableGermlineVariants.add(reportableGermlineVariantWithDriverLikelihood(germlineVariant, 1.0));
                    } else if (filterMinCopyNumberTumor || chordStatus == ChordStatus.HR_DEFICIENT) {
                        reportableGermlineVariants.add(reportableGermlineVariantWithDriverLikelihood(germlineVariant, 0.5));
                    }
                }
            }
        }

        return reportableGermlineVariants;
    }

    @NotNull
    private static List<ReportableGermlineVariant> presentInTumorOnly(@NotNull List<ReportableGermlineVariant> variants) {
        List<ReportableGermlineVariant> variantsInTumor = Lists.newArrayList();
        for (ReportableGermlineVariant variant : variants) {
            if (isPresentInTumor(variant)) {
                variantsInTumor.add(variant);
            }
        }
        return variantsInTumor;
    }

    private static boolean isPresentInTumor(@NotNull ReportableGermlineVariant germlineVariant) {
        return germlineVariant.adjustedCopyNumber() * germlineVariant.adjustedVaf() >= 0.5;
    }

    @NotNull
    private static DriverGermlineVariant reportableGermlineVariantWithDriverLikelihood(
            @NotNull ReportableGermlineVariant germlineVariant, double driverLikelihood) {
        return ImmutableDriverGermlineVariant.builder().variant(germlineVariant).driverLikelihood(driverLikelihood).build();
    }

    @NotNull
    private static GeneCopyNumber lookupGeneCopyNumber(@NotNull List<GeneCopyNumber> allGeneCopyNumbers, @NotNull String gene) {
        for (GeneCopyNumber geneCopyNumber : allGeneCopyNumbers) {
            if (geneCopyNumber.gene().equals(gene)) {
                return geneCopyNumber;
            }
        }

        throw new IllegalStateException("Could not find gene copy number for gene: " + gene);
    }
}
