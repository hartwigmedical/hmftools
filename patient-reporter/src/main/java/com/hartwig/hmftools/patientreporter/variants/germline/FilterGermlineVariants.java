package com.hartwig.hmftools.patientreporter.variants.germline;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordAnalyzer;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.patientreporter.variants.driver.DriverGeneView;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientreporter.variants.ImmutableReportableGermlineVariant;
import com.hartwig.hmftools.patientreporter.variants.ReportableGermlineVariant;

import org.jetbrains.annotations.NotNull;

public final class FilterGermlineVariants {

    private FilterGermlineVariants() {
    }

    @NotNull
    public static List<ReportableGermlineVariant> filterGermlineVariantsForReporting(@NotNull List<GermlineVariant> germlineVariants,
            @NotNull DriverGeneView driverGeneView, @NotNull GermlineReportingModel germlineReportingModel,
            @NotNull List<GeneCopyNumber> allGeneCopyNumbers, @NotNull List<SomaticVariant> variantsToReport,
            @NotNull ChordAnalyzer chordAnalyzer) {
        List<ReportableGermlineVariant> reportableGermlineVariants = Lists.newArrayList();

        Set<String> reportableGermlineGenes = germlineReportingModel.reportableGermlineGenes();
        for (GermlineVariant germlineVariant : germlineVariants) {
            assert germlineVariant.passFilter();

            if (reportableGermlineGenes.contains(germlineVariant.gene())) {
                // Note: Reportable germline genes may not necessarily be present in driverGeneView!
                if (driverGeneView.category(germlineVariant.gene()) == DriverCategory.ONCO) {
                    // Report all germline variants on reportable oncogenes.
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
                    for (SomaticVariant variant : variantsToReport) {
                        if (variant.gene().equals(germlineVariant.gene())) {
                            filterSomaticVariantInSameGene = true;
                        }
                    }

                    boolean filterGermlineVariantInSameGene = false;
                    for (GermlineVariant variant : germlineVariants) {
                        if (variant != germlineVariant && variant.gene().equals(germlineVariant.gene())) {
                            filterGermlineVariantInSameGene = true;
                        }
                    }

                    if (filterBiallelic || filterSomaticVariantInSameGene || filterGermlineVariantInSameGene) {
                        reportableGermlineVariants.add(reportableGermlineVariantWithDriverLikelihood(germlineVariant, 1.0));
                    } else if (filterMinCopyNumberTumor || chordAnalyzer.chordAnalysis().predictedResponseValue()) {
                        reportableGermlineVariants.add(reportableGermlineVariantWithDriverLikelihood(germlineVariant, 0.5));
                    }
                }
            }
        }
        return reportableGermlineVariants;
    }

    @NotNull
    private static ReportableGermlineVariant reportableGermlineVariantWithDriverLikelihood(@NotNull GermlineVariant germlineVariant,
            double driverLikelihood) {
        return ImmutableReportableGermlineVariant.builder().variant(germlineVariant).driverLikelihood(driverLikelihood).build();
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
