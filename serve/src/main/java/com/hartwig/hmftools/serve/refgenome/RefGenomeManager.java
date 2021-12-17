package com.hartwig.hmftools.serve.refgenome;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ImmutableExtractionResult;
import com.hartwig.hmftools.serve.refgenome.liftover.LiftOverAlgo;
import com.hartwig.hmftools.serve.refgenome.liftover.UCSCLiftOver;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RefGenomeManager {

    private static final Logger LOGGER = LogManager.getLogger(RefGenomeManager.class);

    @NotNull
    private final Map<RefGenomeVersion, RefGenomeResource> refGenomeResourceMap;
    @NotNull
    private final ConversionFilter conversionFilter;

    RefGenomeManager(@NotNull final Map<RefGenomeVersion, RefGenomeResource> refGenomeResourceMap) {
        this.refGenomeResourceMap = refGenomeResourceMap;
        this.conversionFilter = new ConversionFilter();
    }

    @NotNull
    public RefGenomeResource pickResourceForKnowledgebase(@NotNull Knowledgebase knowledgebase) {
        return checkedRetrieve(knowledgebase.refGenomeVersion());
    }

    @NotNull
    public IndexedFastaSequenceFile refSequenceForRefGenome(@NotNull RefGenomeVersion version) {
        return checkedRetrieve(version).refSequence();
    }

    public void evaluate() {
        evaluateProteinResolving();

        conversionFilter.reportUnusedFilterEntries();
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

        RefGenomeResource sourceResource = checkedRetrieve(sourceVersion);
        IndexedFastaSequenceFile targetSequence = refSequenceForRefGenome(targetVersion);
        String chainFromSourceToTarget = sourceResource.chainToOtherRefGenomeMap().get(targetVersion);

        LiftOverAlgo liftOverAlgo = UCSCLiftOver.fromChainFile(chainFromSourceToTarget, targetVersion);
        RefGenomeConverter converter = new RefGenomeConverter(sourceVersion, targetVersion, targetSequence, liftOverAlgo);
        ExtractionResult filteredExtraction = conversionFilter.filter(extraction);

        return ImmutableExtractionResult.builder()
                .refGenomeVersion(targetVersion)
                .knownHotspots(converter.convertKnownHotspots(filteredExtraction.knownHotspots()))
                .knownCodons(converter.convertKnownCodons(filteredExtraction.knownCodons()))
                .knownExons(converter.convertKnownExons(filteredExtraction.knownExons()))
                .knownCopyNumbers(converter.convertKnownCopyNumbers(filteredExtraction.knownCopyNumbers()))
                .knownFusionPairs(converter.convertKnownFusionPairs(filteredExtraction.knownFusionPairs()))
                .actionableHotspots(converter.convertActionableHotspots(filteredExtraction.actionableHotspots()))
                .actionableRanges(converter.convertActionableRanges(filteredExtraction.actionableRanges()))
                .actionableGenes(converter.convertActionableGenes(filteredExtraction.actionableGenes()))
                .actionableFusions(converter.convertActionableFusions(filteredExtraction.actionableFusions()))
                .actionableCharacteristics(filteredExtraction.actionableCharacteristics())
                .build();
    }

    @NotNull
    private RefGenomeResource checkedRetrieve(@NotNull RefGenomeVersion version) {
        RefGenomeResource resource = refGenomeResourceMap.get(version);
        if (resource == null) {
            throw new IllegalStateException("No ref genome resources found for ref genome version " + version);
        }
        return resource;
    }

    private void evaluateProteinResolving() {
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
