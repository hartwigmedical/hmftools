package com.hartwig.hmftools.serve.refgenome;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.GeneNameMapping;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ImmutableExtractionResult;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RefGenomeManager {

    private static final Logger LOGGER = LogManager.getLogger(RefGenomeManager.class);

    @NotNull
    private final Map<RefGenomeVersion, RefGenomeResource> refGenomeResourceMap;
    @NotNull
    private final GeneNameMapping geneNameMapping = new GeneNameMapping();

    RefGenomeManager(@NotNull final Map<RefGenomeVersion, RefGenomeResource> refGenomeResourceMap) {
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

    @NotNull
    public Map<RefGenomeVersion, List<ExtractionResult>> makeVersioned(@NotNull List<ExtractionResult> extractions) {
        Map<RefGenomeVersion, List<ExtractionResult>> versionedExtractionMap = Maps.newHashMap();

        for (RefGenomeVersion version : refGenomeResourceMap.keySet()) {
            LOGGER.info("Creating extraction results for ref genome version {}", version);
            List<ExtractionResult> targetExtractions = Lists.newArrayList();
            for (ExtractionResult extraction : extractions) {
                targetExtractions.add(convert(extraction, version));
            }
            versionedExtractionMap.put(version, targetExtractions);
        }

        return versionedExtractionMap;
    }

    @NotNull
    private ExtractionResult convert(@NotNull ExtractionResult extraction, @NotNull RefGenomeVersion targetVersion) {
        RefGenomeVersion sourceVersion = extraction.refGenomeVersion();
        if (sourceVersion == targetVersion) {
            return extraction;
        }

        RefGenomeResource sourceResource = refGenomeResourceMap.get(sourceVersion);
        IndexedFastaSequenceFile targetSequence = refGenomeResourceMap.get(targetVersion).refSequence();

        LiftOver liftOverToTarget = new LiftOver(new File(sourceResource.chainToOtherRefGenomeMap().get(targetVersion)));
        RefGenomeConverter converter =
                new RefGenomeConverter(sourceVersion, targetVersion, targetSequence, liftOverToTarget, geneNameMapping);

        return ImmutableExtractionResult.builder()
                .refGenomeVersion(targetVersion)
                .knownHotspots(converter.convertKnownHotspots(extraction.knownHotspots()))
                .knownCodons(converter.convertKnownCodons(extraction.knownCodons()))
                .knownExons(converter.convertKnownExons(extraction.knownExons()))
                .knownCopyNumbers(converter.convertKnownCopyNumbers(extraction.knownCopyNumbers()))
                .knownFusionPairs(converter.convertKnownFusionPairs(extraction.knownFusionPairs()))
                .actionableHotspots(converter.convertActionableHotspots(extraction.actionableHotspots()))
                .actionableRanges(converter.convertActionableRanges(extraction.actionableRanges()))
                .actionableGenes(converter.convertActionableGenes(extraction.actionableGenes()))
                .actionableFusions(converter.convertActionableFusion(extraction.actionableFusions()))
                .actionableCharacteristics(extraction.actionableCharacteristics())
                .build();
    }
}
