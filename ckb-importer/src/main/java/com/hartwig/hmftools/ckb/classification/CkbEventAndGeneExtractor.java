package com.hartwig.hmftools.ckb.classification;

import com.hartwig.hmftools.ckb.datamodel.variant.Variant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class CkbEventAndGeneExtractor {

    private static final Logger LOGGER = LogManager.getLogger(CkbEventAndGeneExtractor.class);

    private CkbEventAndGeneExtractor() {
    }

    @NotNull
    public static String extractGene(@NotNull Variant variant) {
        String primaryGene = variant.gene().geneSymbol();
        if (primaryGene.equals(CkbConstants.NO_GENE)) {
            return CkbConstants.NO_GENE;
        } else if (CkbConstants.UNMAPPABLE_GENES.contains(primaryGene)) {
            LOGGER.debug("Skipping gene curation for '{}' since gene is unmappable", variant.gene().geneSymbol());
            return primaryGene;
        }

        // TODO Potentially map gene names against ensembl data cache and switch to synonym instead.

        return primaryGene;
    }

    @NotNull
    public static String extractEvent(@NotNull Variant variant) {
        if (isPromiscuousFusion(variant)) {
            return "fusion promiscuous";
        } else if (isFusion(variant)) {
            return removeAllSpaces(variant.variant());
        } else if (variant.variant().contains("exon") && !variant.variant().contains("exon ")) {
            // Some exon variants do not contain a space after the "exon" keyword
            return variant.variant().replace("exon", "exon ");
        } else {
            return variant.variant();
        }
    }

    private static boolean isFusion(@NotNull Variant variant) {
        return variant.impact() != null && variant.impact().equals("fusion");
    }

    private static boolean isPromiscuousFusion(@NotNull Variant variant) {
        return isFusion(variant) && variant.variant().equals("fusion");
    }

    @NotNull
    private static String removeAllSpaces(@NotNull String value) {
        return value.replaceAll("\\s+", "");
    }
}
