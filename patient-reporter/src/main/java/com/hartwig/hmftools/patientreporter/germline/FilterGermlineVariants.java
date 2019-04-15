package com.hartwig.hmftools.patientreporter.germline;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class FilterGermlineVariants {
    private static final Logger LOGGER = LogManager.getLogger(FilterGermlineVariants.class);

    private FilterGermlineVariants() {
    }

    @NotNull
    public static List<GermlineVariant> filteringReportedGermlineVariant(List<GermlineVariant> germlineVariants,
            @NotNull GermlineGenesReporting germlineGenesReporting, @NotNull GeneModel panelGeneModel,
            @NotNull List<GeneCopyNumber> geneCopyNumbers) {
        List<GermlineVariant> filteredGermlineVariant = Lists.newArrayList();
        Set<String> reportingGenes = germlineGenesReporting.germlineGenes();
        Set<String> notifyGenes = germlineGenesReporting.germlineGenesNotify();

        for (GermlineVariant germlineVariant : germlineVariants) {
            if (germlineVariant.passFilter()) {
                if (panelGeneModel.geneDriverCategoryMap().get(germlineVariant.gene()) == DriverCategory.ONCO) { // use all genes
                    filteredGermlineVariant.add(germlineVariant);
                } else if (panelGeneModel.geneDriverCategoryMap().get(germlineVariant.gene()) == DriverCategory.TSG) { // filter genes
                    filteredGermlineVariant.add(germlineVariant);


                    //                    if (reportingGenes.contains(germlineVariant.gene())) {
                    //                        LOGGER.info("filtering");
                    //
                    //                        if (germlineVariant.biallelic()) { // filter for biallelic
                    //                            for (GeneCopyNumber geneCopyNumber: geneCopyNumbers) {
                    //                                if (Doubles.equal(geneCopyNumber.maxCopyNumber(), 1)) { // filter for gene copy number
                    //                                    filteredGermlineVariant.add(germlineVariant);
                    //                                }
                    //                            }
                    //                        }
                    //                    }
                }
            }
        }

        return filteredGermlineVariant;
    }
}
