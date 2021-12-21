package com.hartwig.hmftools.serve.transvar;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.GeneNameMapping;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class TransvarCurator {

    private static final Logger LOGGER = LogManager.getLogger(TransvarCurator.class);

    private static final Map<String, String> PROTEIN_ANNOTATION_MAPPING = Maps.newHashMap();
    private static final Set<String> GENES_FOR_WHICH_TO_SKIP_MAPPING = Sets.newHashSet();

    static {
        // Transvar can't interpret start-lost so we map it to another arbitrary mutation.
        PROTEIN_ANNOTATION_MAPPING.put("M1?", "M1I");

        // These genes work fine in transvar and should not be mapped to our old HMF v37 name.
        GENES_FOR_WHICH_TO_SKIP_MAPPING.add("CCDC186");
    }

    @NotNull
    private final GeneNameMapping geneNameMapping = new GeneNameMapping();

    public TransvarCurator() {
    }

    @NotNull
    public String curateGene(@NotNull String gene) {
        /// Transvar uses its own gene model which is different from HGNC and also slightly from our old v37 model.
        String transvarGene = gene;
        if (GENES_FOR_WHICH_TO_SKIP_MAPPING.contains(gene)) {
            transvarGene = gene;
            LOGGER.debug("Skipping mapping gene '{}' for transvar", gene);
        } else if (geneNameMapping.hasNewGene(gene)) {
            transvarGene = geneNameMapping.getOldName(gene);
        }

        return transvarGene;
    }

    @NotNull
    public String curateProteinAnnotation(@NotNull String proteinAnnotation) {
        String modifiedProteinAnnotation = proteinAnnotation;
        if (PROTEIN_ANNOTATION_MAPPING.containsKey(proteinAnnotation)) {
            modifiedProteinAnnotation = PROTEIN_ANNOTATION_MAPPING.get(proteinAnnotation);
            LOGGER.debug("Mapping protein annotation '{}' to '{}'", proteinAnnotation, modifiedProteinAnnotation);
        }
        return modifiedProteinAnnotation;
    }
}
