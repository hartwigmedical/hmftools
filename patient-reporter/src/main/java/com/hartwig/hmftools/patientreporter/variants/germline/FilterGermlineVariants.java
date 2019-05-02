package com.hartwig.hmftools.patientreporter.variants.germline;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;

import org.jetbrains.annotations.NotNull;

public final class FilterGermlineVariants {

    private FilterGermlineVariants() {
    }

    @NotNull
    public static List<GermlineVariant> filterGermlineVariantsForReporting(List<GermlineVariant> germlineVariants,
            @NotNull GermlineReportingModel germlineReportingModel, @NotNull Map<String, DriverCategory> driverCategoryPerGeneMap,
            @NotNull List<GeneCopyNumber> allGeneCopyNumbers, @NotNull List<EnrichedSomaticVariant> variantsToReport) {
        List<GermlineVariant> filteredGermlineVariant = Lists.newArrayList();

        Set<String> reportingGermlineGenes = germlineReportingModel.reportableGermlineGenes();
        for (GermlineVariant germlineVariant : germlineVariants) {
            assert germlineVariant.passFilter();

            if (reportingGermlineGenes.contains(germlineVariant.gene())) {
                if (driverCategoryPerGeneMap.get(germlineVariant.gene()) == DriverCategory.ONCO) {
                    // Report all germline variants on reportable oncogenes.
                    filteredGermlineVariant.add(germlineVariant);
                } else {
                    // Only report germline variants on TSGs if there is a 2nd hit.
                    boolean filterBiallelic = germlineVariant.biallelic();

                    boolean filterMinCopyNumberTumor = false;
                    GeneCopyNumber geneCopyNumber = lookupGeneCopyNumber(allGeneCopyNumbers, germlineVariant.gene());
                    if (Math.round(geneCopyNumber.minCopyNumber()) < 2) {
                        filterMinCopyNumberTumor = true;
                    }

                    boolean filterSomaticVariantInSameGene = false;
                    for (EnrichedSomaticVariant variant : variantsToReport) {
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
