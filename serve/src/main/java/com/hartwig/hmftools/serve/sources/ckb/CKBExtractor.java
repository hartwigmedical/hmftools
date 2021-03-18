package com.hartwig.hmftools.serve.sources.ckb;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.ckb.classification.EventTypeExtractor;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.actionability.signature.ActionableSignature;
import com.hartwig.hmftools.serve.extraction.EventExtractor;
import com.hartwig.hmftools.serve.extraction.ExtractionFunctions;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ImmutableExtractionResult;
import com.hartwig.hmftools.serve.sources.vicc.ViccUtil;
import com.hartwig.hmftools.serve.util.ProgressTracker;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CKBExtractor {

    private static final Logger LOGGER = LogManager.getLogger(CKBExtractor.class);

    @NotNull
    private final EventExtractor eventExtractor;
    @NotNull
    private final ActionableEvidenceFactory actionableEvidenceFactory;

    public CKBExtractor(@NotNull final EventExtractor eventExtractor, @NotNull final ActionableEvidenceFactory actionableEvidenceFactory) {
        this.eventExtractor = eventExtractor;
        this.actionableEvidenceFactory = actionableEvidenceFactory;
    }

    @NotNull
    public ExtractionResult extract(@NotNull List<CkbEntry> ckbEntries) {
        Map<CkbEntry, ExtractResult> resultsPerEntry = Maps.newHashMap();

        ProgressTracker tracker = new ProgressTracker("CKB", ckbEntries.size());
        for (CkbEntry entry : ckbEntries) {
            //   resultsPerEntry.put(entry, extractSingleEntry(entry));
            LOGGER.info(entry.profileName());
            LOGGER.info(entry.variants());

                for (Variant variant : entry.variants()) {
                    LOGGER.info("profileName: {} ", entry.profileName());
                    LOGGER.info("eventType: {}", EventTypeExtractor.classify(variant.gene().geneSymbol(), variant.variant()));
                    eventExtractor.extract(variant.gene().geneSymbol(),
                            variant.gene().canonicalTranscript(),
                            EventTypeExtractor.classify(variant.gene().geneSymbol(), variant.variant()),
                            variant.impact());

                }


            tracker.update();
        }

        // actionableEvidenceFactory.evaluateCuration();

        CKBUtils.printExtractionResults(resultsPerEntry);


        ImmutableExtractionResult.Builder outputBuilder = ImmutableExtractionResult.builder()
                .knownHotspots(Sets.newHashSet())
                .knownCodons(Sets.newHashSet())
                .knownExons(Sets.newHashSet())
                .knownCopyNumbers(Sets.newHashSet())
                .knownFusionPairs(Sets.newHashSet());

        addActionability(outputBuilder);

        return ExtractionFunctions.consolidateActionableEvents(outputBuilder.build());
    }

    private static void addActionability(@NotNull ImmutableExtractionResult.Builder outputBuilder) {
        Set<ActionableHotspot> actionableHotspots = Sets.newHashSet();
        Set<ActionableRange> actionableRanges = Sets.newHashSet();
        Set<ActionableGene> actionableGenes = Sets.newHashSet();
        Set<ActionableFusion> actionableFusions = Sets.newHashSet();
        Set<ActionableSignature> actionableSignatures = Sets.newHashSet();

        outputBuilder.actionableHotspots(actionableHotspots);
        outputBuilder.actionableRanges(actionableRanges);
        outputBuilder.actionableGenes(actionableGenes);
        outputBuilder.actionableFusions(actionableFusions);
        outputBuilder.actionableSignatures(actionableSignatures);
    }

}
