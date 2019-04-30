package com.hartwig.hmftools.patientreporter.germline;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingChoice;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;

import org.jetbrains.annotations.NotNull;

public final class FilterGermlineVariants {

    private FilterGermlineVariants() {
    }

    @NotNull
    public static List<GermlineVariant> filterGermlineVariantsForReporting(List<GermlineVariant> germlineVariants,
            @NotNull Map<String, Boolean> germlineGenesReporting, @NotNull Map<String, DriverCategory> driverCategoryPerGeneMap,
            @NotNull List<GeneCopyNumber> allGeneCopyNumbers,
            @NotNull List<EnrichedSomaticVariant> variantsToReport, @NotNull LimsGermlineReportingChoice choiceGermlineReporting) {
        List<GermlineVariant> filteredGermlineVariant = Lists.newArrayList();

        Set<String> reportingGermlineGenes = germlineGenesReporting.keySet();

        for (GermlineVariant germlineVariant : germlineVariants) {
            boolean filterBiallelic = false;
            boolean filterMinCopyNumberTumor = false;
            boolean filterSomaticVariantInSameGene = false;
            if (germlineVariant.passFilter() && reportingGermlineGenes.contains(germlineVariant.gene()) && !choiceGermlineReporting.equals(
                    LimsGermlineReportingChoice.UNKNOWN)) {
                if (driverCategoryPerGeneMap.get(germlineVariant.gene()) == DriverCategory.ONCO) { // use all genes
                    filteredGermlineVariant.add(germlineVariant);
                } else { // filter genes
                    if (germlineVariant.biallelic()) { // variant is biallelic (2nd hit CNV)
                        filterBiallelic = true;
                    }
                    for (GeneCopyNumber geneCopyNumber : allGeneCopyNumbers) { // min copy number in tumer = 1 (2nd hit SV)
                        if (geneCopyNumber.gene().equals(germlineVariant.gene())) { // filter for gene copy number
                            if (Math.round(geneCopyNumber.minCopyNumber()) == 1) {
                                filterMinCopyNumberTumor = true;
                            }
                        }
                    }
                    for (EnrichedSomaticVariant variant : variantsToReport) { // filter for gene copy number
                        if (variant.gene().equals(germlineVariant.gene())) {
                            filterSomaticVariantInSameGene = true;
                        }
                    }
                    if (filterBiallelic || filterMinCopyNumberTumor || filterSomaticVariantInSameGene) {
                        filteredGermlineVariant.add(germlineVariant);
                    }
                }
            }
        }
        return filteredGermlineVariant;
    }
}
