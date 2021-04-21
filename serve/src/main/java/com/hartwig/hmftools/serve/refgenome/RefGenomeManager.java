package com.hartwig.hmftools.serve.refgenome;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ImmutableExtractionResult;

import org.apache.commons.compress.utils.Lists;
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

    @NotNull
    public Map<RefGenomeVersion, List<ExtractionResult>> makeVersioned(@NotNull List<ExtractionResult> extractions) {
        Map<RefGenomeVersion, List<ExtractionResult>> versionedExtractionMap = Maps.newHashMap();

        for (RefGenomeVersion version : refGenomeResourceMap.keySet()) {
            LOGGER.info("Creating extraction results for ref genome version {}", version);
            List<ExtractionResult> targetExtractions = Lists.newArrayList();
            for (ExtractionResult extraction : extractions) {
//                if (extraction.refGenomeVersion() == version) {
//                    targetExtractions.add(extraction);
//                } else {
                    targetExtractions.add(convert(extraction, version));
//                }
            }
            versionedExtractionMap.put(version, targetExtractions);
        }

        return versionedExtractionMap;
    }

    @NotNull
    private ExtractionResult convert(@NotNull ExtractionResult extraction, @NotNull RefGenomeVersion targetVersion) {
        // TODO: Convert hotspots: position, chromosome, gene. Check if ref is identical.
        // TODO: Convert codons: gene, chromosome, start, end.
        // TODO: Convert exons: gene, chromosome, start, end
        // TODO: Convert known copy numbers: gene
        // TODO: Convert known fusions: geneUp, geneDown
        // TODO: Actionable follows known

        return ImmutableExtractionResult.builder().from(extraction).refGenomeVersion(RefGenomeVersion.V37).build();
    }
}
