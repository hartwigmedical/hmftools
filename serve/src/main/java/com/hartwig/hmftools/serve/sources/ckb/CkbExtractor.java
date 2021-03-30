package com.hartwig.hmftools.serve.sources.ckb;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.extraction.ActionableEventFactory;
import com.hartwig.hmftools.serve.extraction.EventExtractor;
import com.hartwig.hmftools.serve.extraction.EventExtractorOutput;
import com.hartwig.hmftools.serve.extraction.ExtractionFunctions;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ImmutableExtractionResult;
import com.hartwig.hmftools.serve.sources.iclusion.ActionableTrial;
import com.hartwig.hmftools.serve.util.ProgressTracker;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CkbExtractor {

    private static final Logger LOGGER = LogManager.getLogger(CkbExtractor.class);

    @NotNull
    private final EventExtractor eventExtractor;
    @NotNull
    private final ActionableEvidenceFactory actionableEvidenceFactory;

    public CkbExtractor(@NotNull final EventExtractor eventExtractor, @NotNull final ActionableEvidenceFactory actionableEvidenceFactory) {
        this.eventExtractor = eventExtractor;
        this.actionableEvidenceFactory = actionableEvidenceFactory;
    }

    @NotNull
    public ExtractionResult extract(@NotNull List<CkbEntry> ckbEntries, @NotNull EventType type) {
        List<ExtractionResult> extractions = Lists.newArrayList();
        Set<ActionableEvent> actionableEvents = Sets.newHashSet();
        List<EventExtractorOutput> eventExtractions = Lists.newArrayList();

        ProgressTracker tracker = new ProgressTracker("CKB", ckbEntries.size());
        for (CkbEntry entry : ckbEntries) {

            if (entry.variants().size() == 1) {
                eventExtractions.add(eventExtractor.extract(entry.variants().get(0).gene().geneSymbol(),
                        entry.variants().get(0).gene().canonicalTranscript(),
                        type,
                        entry.variants().get(0).variant()));
                actionableEvents = actionableEvidenceFactory.toActionableEvents(entry);

                if (type == EventType.UNKNOWN) {
                    LOGGER.warn("No event type known for '{}' on '{}'",
                            entry.variants().get(0).variant(),
                            entry.variants().get(0).gene().geneSymbol());
                }

            }

            extractions.add(toExtractionResult(actionableEvents, eventExtractions));

            tracker.update();
        }

        // actionableEvidenceFactory.evaluateCuration();

        // CkbUtils.printExtractionResults(resultsPerEntry);

        ImmutableExtractionResult.Builder outputBuilder = ImmutableExtractionResult.builder()
                .knownHotspots(Sets.newHashSet())
                .knownCodons(Sets.newHashSet())
                .knownExons(Sets.newHashSet())
                .knownCopyNumbers(Sets.newHashSet())
                .knownFusionPairs(Sets.newHashSet());

        return ExtractionFunctions.merge(extractions);

    }

    @NotNull
    private static ExtractionResult toExtractionResult(@NotNull Set<ActionableEvent> actionableTrials,
            @NotNull List<EventExtractorOutput> eventExtractions) {
        Set<ActionableHotspot> actionableHotspots = Sets.newHashSet();
        Set<ActionableRange> actionableRanges = Sets.newHashSet();
        Set<ActionableGene> actionableGenes = Sets.newHashSet();
        Set<ActionableFusion> actionableFusions = Sets.newHashSet();
        Set<ActionableCharacteristic> actionableCharacteristics = Sets.newHashSet();

        for (ActionableEvent trial : actionableTrials) {
            for (EventExtractorOutput extraction : eventExtractions) {
                actionableHotspots.addAll(ActionableEventFactory.toActionableHotspots(trial, extraction.hotspots()));
                actionableRanges.addAll(ActionableEventFactory.toActionableRanges(trial, extraction.codons()));
                actionableRanges.addAll(ActionableEventFactory.toActionableRanges(trial, extraction.exons()));

                if (extraction.geneLevelEvent() != null) {
                    actionableGenes.add(ActionableEventFactory.geneLevelEventToActionableGene(trial, extraction.geneLevelEvent()));
                }

                if (extraction.knownCopyNumber() != null) {
                    actionableGenes.add(ActionableEventFactory.copyNumberToActionableGene(trial, extraction.knownCopyNumber()));
                }

                if (extraction.knownFusionPair() != null) {
                    actionableFusions.add(ActionableEventFactory.toActionableFusion(trial, extraction.knownFusionPair()));
                }

                if (extraction.characteristic() != null) {
                    actionableCharacteristics.add(ActionableEventFactory.toActionableCharacteristic(trial, extraction.characteristic()));
                }
            }
        }

        return ImmutableExtractionResult.builder()
                .actionableHotspots(actionableHotspots)
                .actionableRanges(actionableRanges)
                .actionableGenes(actionableGenes)
                .actionableFusions(actionableFusions)
                .actionableCharacteristics(actionableCharacteristics)
                .build();
    }

}
