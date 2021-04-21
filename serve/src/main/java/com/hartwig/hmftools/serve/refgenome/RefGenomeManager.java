package com.hartwig.hmftools.serve.refgenome;

import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class RefGenomeManager {

    private static final Logger LOGGER = LogManager.getLogger(RefGenomeManager.class);

    @NotNull
    private final Map<RefGenomeVersion, RefGenomeResource> refGenomeResourceMap;

    RefGenomeManager(@NotNull final Map<RefGenomeVersion, RefGenomeResource> refGenomeResourceMap) {
        this.refGenomeResourceMap = refGenomeResourceMap;
    }

    @NotNull
    public Set<RefGenomeVersion> versions() {
        return refGenomeResourceMap.keySet();
    }

    @NotNull
    public RefGenomeResource pickResourceForVersion(@NotNull RefGenomeVersion version) {
        RefGenomeResource resource = refGenomeResourceMap.get(version);
        if (resource == null) {
            throw new IllegalStateException("No ref genome resource found for version " + version);
        }
        return resource;
    }

    @NotNull
    public RefGenomeResource pickResourceForKnowledgebase(@NotNull Knowledgebase knowledgebase) {
        RefGenomeResource resource = refGenomeResourceMap.get(knowledgebase.refGenomeVersion());
        if (resource == null) {
            throw new IllegalStateException("No ref genome resources found for knowledgebase " + knowledgebase + " with version "
                    + knowledgebase.refGenomeVersion());
        }
        return resource;
    }

    public void evaluateProteinResolving() {
        for (Map.Entry<RefGenomeVersion, RefGenomeResource> entry : refGenomeResourceMap.entrySet()) {
            RefGenomeVersion version = entry.getKey();
            RefGenomeResource resource = entry.getValue();
            Set<String> unresolvedProteinAnnotations = resource.proteinResolver().unresolvedProteinAnnotations();
            if (!unresolvedProteinAnnotations.isEmpty()) {
                LOGGER.warn("Protein resolver {} could not resolve {} protein annotations", version, unresolvedProteinAnnotations.size());
                for (String unresolvedProteinAnnotation : unresolvedProteinAnnotations) {
                    LOGGER.warn("Protein resolver {} could not resolve protein annotation '{}'", version, unresolvedProteinAnnotation);
                }
            } else {
                LOGGER.debug("Protein resolver {} observed no issues when resolving hotspots", version);
            }
        }
    }
}
