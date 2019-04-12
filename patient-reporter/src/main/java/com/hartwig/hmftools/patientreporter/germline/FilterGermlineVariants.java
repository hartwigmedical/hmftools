package com.hartwig.hmftools.patientreporter.germline;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class FilterGermlineVariants {
    private static final Logger LOGGER = LogManager.getLogger(FilterGermlineVariants.class);

    private FilterGermlineVariants() {
    }

    @NotNull
    public static List<GermlineVariant> filteringReportedGermlineVariant(List<GermlineVariant> germlineVariants,
            @NotNull GermlineGenesReporting germlineGenesReporting) {
        List<GermlineVariant> filteredGermlineVariant = Lists.newArrayList();
        Set<String> reportingGenes = germlineGenesReporting.germlineGenes();
        Set<String> notifyGenes = germlineGenesReporting.germlineGenesNotify();

        for (GermlineVariant germlineVariant : germlineVariants) {
            if (reportingGenes.contains(germlineVariant.gene())) {
                if (germlineVariant.biallelic()) {
                    filteredGermlineVariant.add(germlineVariant);
                }
            }
        }

        LOGGER.info("filtered germline variants: ");
        LOGGER.info(filteredGermlineVariant);
        return filteredGermlineVariant;
    }
}
