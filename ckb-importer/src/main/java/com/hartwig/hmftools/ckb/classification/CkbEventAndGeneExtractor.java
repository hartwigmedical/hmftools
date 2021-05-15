package com.hartwig.hmftools.ckb.classification;

import com.hartwig.hmftools.ckb.datamodel.variant.Variant;
import com.hartwig.hmftools.common.genome.refgenome.GeneNameMapping;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CkbEventAndGeneExtractor {

    private static final Logger LOGGER = LogManager.getLogger(CkbEventAndGeneExtractor.class);

    @NotNull
    private final GeneNameMapping geneNameMapping;

    public CkbEventAndGeneExtractor() {
        this.geneNameMapping = GeneNameMapping.loadFromEmbeddedResource();
    }

    @NotNull
    public String extractGene(@NotNull Variant variant) {
        if (isFusion(variant) && !isPromiscuousFusion(variant)) {
            return variant.fullName();
        } else if (variant.gene().geneSymbol().equals(CkbConstants.NO_GENE)) {
            return CkbConstants.NO_GENE;
        } else if (CkbConstants.UNMAPPABLE_GENES.contains(variant.gene().geneSymbol())) {
            LOGGER.debug("Skipping gene curation for '{}' since gene is unmappable", variant.gene().geneSymbol());
            return variant.gene().geneSymbol();
        } else {
            String primaryGene = variant.gene().geneSymbol();
            if (!geneNameMapping.isValidV38Gene(primaryGene)) {
                for (String synonym : variant.gene().synonyms()) {
                    if (geneNameMapping.isValidV38Gene(synonym)) {
                        LOGGER.debug("Swapping CKB gene '{}' with synonym '{}'", primaryGene, synonym);
                        return synonym;
                    }
                }

                // Only warn in case the primary gene is a real gene in the first place.
                if (!CkbConstants.NON_EXISTING_GENES.contains(primaryGene)) {
                    LOGGER.warn("Could not find synonym for '{}' that exists in HMF v38 gene model", primaryGene);
                }
            }
            return primaryGene;
        }
    }

    @NotNull
    public String extractEvent(@NotNull Variant variant) {
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
