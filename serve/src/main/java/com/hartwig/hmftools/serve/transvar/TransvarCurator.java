package com.hartwig.hmftools.serve.transvar;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.GeneNameMapping;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class TransvarCurator {

    private static final Logger LOGGER = LogManager.getLogger(TransvarCurator.class);
    private static final Map<String, String> PROTEIN_ANNOTATION_MAPPING = Maps.newHashMap();

    static {
        // Transvar can't interpret start-lost so we map it to another mutation.
        PROTEIN_ANNOTATION_MAPPING.put("M1?", "M1I");
    }

    @NotNull
    private final RefGenomeVersion refGenomeVersion;
    @NotNull
    private final GeneNameMapping geneNameMapping = GeneNameMapping.loadFromEmbeddedResource();

    public TransvarCurator(@NotNull final RefGenomeVersion refGenomeVersion) {
        this.refGenomeVersion = refGenomeVersion;
    }

    @NotNull
    public String curateGene(@NotNull String gene) {
        // Somehow can't manage to configure transvar to use typical v38 gene names.
        String v37Gene = gene;
        if (refGenomeVersion == RefGenomeVersion.V38) {
            if (geneNameMapping.isValidV38Gene(gene)) {
                v37Gene = geneNameMapping.v37Gene(gene);
            } else {
                LOGGER.warn("Could not map gene '{}' to v37 for transvar!", gene);
            }
        }
        return v37Gene;
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
