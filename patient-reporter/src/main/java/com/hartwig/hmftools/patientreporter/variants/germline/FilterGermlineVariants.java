package com.hartwig.hmftools.patientreporter.variants.germline;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariant;
import com.hartwig.hmftools.patientreporter.variants.ImmutableReportableGermlineVariantExtended;
import com.hartwig.hmftools.patientreporter.variants.ReportableGermlineVariantExtended;

import org.jetbrains.annotations.NotNull;

public final class FilterGermlineVariants {

    private FilterGermlineVariants() {
    }

    @NotNull
    public static List<ReportableGermlineVariantExtended> filterGermlineVariantsForReporting(
            @NotNull List<ReportableGermlineVariant> germlineVariants, @NotNull DriverGenePanel driverGenePanel,
            @NotNull GermlineReportingModel germlineReportingModel, @NotNull List<GeneCopyNumber> allGeneCopyNumbers,
            @NotNull List<SomaticVariant> variantsToReport, @NotNull ChordStatus chordStatus) {
        List<ReportableGermlineVariantExtended> reportableGermlineVariants = Lists.newArrayList();

        Set<String> reportableGermlineGenes = germlineReportingModel.reportableGermlineGenes();
        for (ReportableGermlineVariant germlineVariant : germlineVariants) {
            if (reportableGermlineGenes.contains(germlineVariant.gene())) {
                // Note: Reportable germline genes may not necessarily be present in driverGeneView!
                if (driverGenePanel.category(germlineVariant.gene()) == DriverCategory.ONCO) {
                    // Report all germline variants on reportable oncogenes.
                    if (filterSubclonalGermlineVariants(germlineVariant)) {
                        reportableGermlineVariants.add(reportableGermlineVariantWithDriverLikelihood(germlineVariant, 1.0));
                    }
                } else {
                    // Only report germline variants on TSGs if there is a 2nd hit or CHORD suggests HRD
                    boolean filterBiallelic = germlineVariant.biallelic();

                    boolean filterMinCopyNumberTumor = false;
                    GeneCopyNumber geneCopyNumber = lookupGeneCopyNumber(allGeneCopyNumbers, germlineVariant.gene());
                    if (Math.round(geneCopyNumber.minCopyNumber()) <= 1 && (Math.round(germlineVariant.adjustedCopyNumber()) >= 2)) {
                        filterMinCopyNumberTumor = true;
                    }

                    boolean filterSomaticVariantInSameGene = false;
                    for (SomaticVariant variant : variantsToReport) {
                        if (variant.gene().equals(germlineVariant.gene())) {
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
                        if (filterSubclonalGermlineVariants(germlineVariant)) {
                            reportableGermlineVariants.add(reportableGermlineVariantWithDriverLikelihood(germlineVariant, 1.0));
                        }
                    } else if (filterMinCopyNumberTumor || chordStatus == ChordStatus.HRD) {
                        if (filterSubclonalGermlineVariants(germlineVariant)) {
                            reportableGermlineVariants.add(reportableGermlineVariantWithDriverLikelihood(germlineVariant, 0.5));
                        }
                    }
                }
            }
        }

        return reportableGermlineVariants;
    }

    @NotNull
    private static ReportableGermlineVariantExtended reportableGermlineVariantWithDriverLikelihood(
            @NotNull ReportableGermlineVariant germlineVariant, double driverLikelihood) {
        return ImmutableReportableGermlineVariantExtended.builder().variant(germlineVariant).driverLikelihood(driverLikelihood).build();
    }

    private static boolean filterSubclonalGermlineVariants(@NotNull ReportableGermlineVariant germlineVariant) {
        //Filter germline variants out when it is subclonal
        if (germlineVariant.adjustedVaf() >= 0.5) {
            return true;
        } else {
            return false;
        }
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
