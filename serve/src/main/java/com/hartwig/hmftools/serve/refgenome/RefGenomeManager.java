package com.hartwig.hmftools.serve.refgenome;

import java.io.FileNotFoundException;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.ServeConfig;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolverFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class RefGenomeManager {

    private static final Logger LOGGER = LogManager.getLogger(RefGenomeManager.class);

    @NotNull
    private final Map<RefGenomeVersion, RefGenomeResource> refGenomeResourceMap;

    @NotNull
    public static RefGenomeManager buildFromServeConfig(@NotNull ServeConfig config) throws FileNotFoundException {
        Map<RefGenomeVersion, RefGenomeResource> refGenomeResourceMap = Maps.newHashMap();
        refGenomeResourceMap.put(RefGenomeVersion.V37, buildRefGenomeResource37(config));
        refGenomeResourceMap.put(RefGenomeVersion.V38, buildRefGenomeResource38(config));
        return new RefGenomeManager(refGenomeResourceMap);
    }

    @NotNull
    private static RefGenomeResource buildRefGenomeResource37(@NotNull ServeConfig config) throws FileNotFoundException {
        String fastaFile37 = config.refGenome37FastaFile();
        LOGGER.info("Creating ref genome resource for V37 using fasta {}", fastaFile37);
        Map<String, HmfTranscriptRegion> transcriptMap37 = HmfGenePanelSupplier.allGenesMap37();
        ProteinResolver proteinResolver37 = config.skipHotspotResolving()
                ? ProteinResolverFactory.dummy()
                : ProteinResolverFactory.transvarWithRefGenome(RefGenomeVersion.V37, fastaFile37, transcriptMap37);

        return ImmutableRefGenomeResource.builder()
                .fastaFile(fastaFile37)
                .canonicalTranscriptPerGeneMap(transcriptMap37)
                .putChainToOtherRefGenomeMap(RefGenomeVersion.V38, config.refGenome37To38Chain())
                .proteinResolver(proteinResolver37)
                .build();
    }

    @NotNull
    private static RefGenomeResource buildRefGenomeResource38(@NotNull ServeConfig config) throws FileNotFoundException {
        String fastaFile38 = config.refGenome38FastaFile();
        LOGGER.info("Creating ref genome resource for V38 using fasta {}", fastaFile38);
        Map<String, HmfTranscriptRegion> transcriptMap38 = HmfGenePanelSupplier.allGenesMap38();
        ProteinResolver proteinResolver38 = config.skipHotspotResolving()
                ? ProteinResolverFactory.dummy()
                : ProteinResolverFactory.transvarWithRefGenome(RefGenomeVersion.V38, fastaFile38, transcriptMap38);

        return ImmutableRefGenomeResource.builder()
                .fastaFile(fastaFile38)
                .canonicalTranscriptPerGeneMap(transcriptMap38)
                .putChainToOtherRefGenomeMap(RefGenomeVersion.V37, config.refGenome38To37Chain())
                .proteinResolver(proteinResolver38)
                .build();
    }

    private RefGenomeManager(@NotNull final Map<RefGenomeVersion, RefGenomeResource> refGenomeResourceMap) {
        this.refGenomeResourceMap = refGenomeResourceMap;
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
