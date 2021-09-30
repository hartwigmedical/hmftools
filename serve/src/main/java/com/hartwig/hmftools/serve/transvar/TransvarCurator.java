package com.hartwig.hmftools.serve.transvar;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.genepanel.GeneNameMapping37to38;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class TransvarCurator {

    private static final Logger LOGGER = LogManager.getLogger(TransvarCurator.class);

    private static final Map<String, String> PROTEIN_ANNOTATION_MAPPING = Maps.newHashMap();
    private static final Set<String> GENES_FOR_WHICH_TO_SKIP_38_MAPPING = Sets.newHashSet();
    private static final Map<String, String> MANUAL_GENE_MAPPING_38_TO_37 = Maps.newHashMap();

    static {
        // Transvar can't interpret start-lost so we map it to another arbitrary mutation.
        PROTEIN_ANNOTATION_MAPPING.put("M1?", "M1I");

        // These genes have to be mapped for transvar specifically since the HMF v37 gene name does not exist in transvar.
        MANUAL_GENE_MAPPING_38_TO_37.put("EPOP", "C17orf96");
        MANUAL_GENE_MAPPING_38_TO_37.put("NAPB", "SEPT9");

        // These genes work fine in transvar and should not be mapped to the HMF v37 name.
        GENES_FOR_WHICH_TO_SKIP_38_MAPPING.add("CCDC186");
    }

    @NotNull
    private final RefGenomeVersion refGenomeVersion;
    @NotNull
    private final GeneNameMapping37to38 geneNameMapping = GeneNameMapping37to38.loadFromEmbeddedResource();

    public TransvarCurator(@NotNull final RefGenomeVersion refGenomeVersion) {
        this.refGenomeVersion = refGenomeVersion;
    }

    @NotNull
    public String curateGene(@NotNull String gene) {
        // Somehow can't manage to configure transvar to use typical v38 gene names.
        String v37Gene = gene;
        if (refGenomeVersion == RefGenomeVersion.V38) {
            if (MANUAL_GENE_MAPPING_38_TO_37.containsKey(gene)) {
                v37Gene = MANUAL_GENE_MAPPING_38_TO_37.get(gene);
                LOGGER.debug("Manually mapping gene '{}' to '{}' for transvar", gene, v37Gene);
            } else if (GENES_FOR_WHICH_TO_SKIP_38_MAPPING.contains(gene)) {
                v37Gene = gene;
                LOGGER.debug("Skipping mapping gene '{}' for transvar", gene);
            } else if (geneNameMapping.isValidV38Gene(gene)) {
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
