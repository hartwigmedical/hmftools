package com.hartwig.hmftools.serve.sources.ckb;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.common.refseq.RefSeq;
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
import com.hartwig.hmftools.ckb.classification.EventExtractorCuration;
import com.hartwig.hmftools.serve.util.ProgressTracker;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

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

    @Nullable
    public String extractCanonicalTranscript(@NotNull String refseqToMatch, @NotNull List<RefSeq> refSeqMatchFile) {
        for (RefSeq refSeq : refSeqMatchFile) {
            if (refSeq.dbPrimaryAcc().equals(refseqToMatch)) {
                return refSeq.transcriptId();
            }
        }
        return null;
    }

    @NotNull
    public ExtractionResult extract(@NotNull List<CkbEntry> ckbEntries, @NotNull List<RefSeq> refSeqMatchFile) {
        List<ExtractionResult> extractions = Lists.newArrayList();

        ProgressTracker tracker = new ProgressTracker("CKB", ckbEntries.size());
        for (CkbEntry entry : ckbEntries) {
            if (entry.variants().size() == 1) {
                String geneSymbol = EventExtractorCuration.extractGene(entry);
                String event = EventExtractorCuration.extractEvent(entry);

                EventExtractorOutput eventExtractorOutput = eventExtractor.extract(geneSymbol,
                        extractCanonicalTranscript(entry.variants().get(0).gene().canonicalTranscript(), refSeqMatchFile),
                        entry.type(),
                        event);
                Set<ActionableEvent> actionableEvents = actionableEvidenceFactory.toActionableEvents(entry);

                CkbExtractorResult ckbExtractorResult = toResult(eventExtractorOutput, actionableEvents);

                extractions.add(toExtractionResult(actionableEvents, ckbExtractorResult));

                if (entry.type() == EventType.UNKNOWN) {
                    LOGGER.warn("No event type known for '{}' on '{}'",
                            entry.variants().get(0).variant(),
                            entry.variants().get(0).gene().geneSymbol());
                }
            }
            tracker.update();
        }

        actionableEvidenceFactory.evaluateCuration();

        return ExtractionFunctions.merge(extractions);
    }

    @NotNull
    private CkbExtractorResult toResult(@NotNull EventExtractorOutput eventExtractorOutput,
            @NotNull Set<ActionableEvent> actionableEvents) {

        return ImmutableCkbExtractorResult.builder()
                .hotspots(eventExtractorOutput.hotspots())
                .codons(eventExtractorOutput.codons())
                .exons(eventExtractorOutput.exons())
                .geneLevelEvent(eventExtractorOutput.geneLevelEvent())
                .knownCopyNumber(eventExtractorOutput.knownCopyNumber())
                .knownFusionPair(eventExtractorOutput.knownFusionPair())
                .characteristic(eventExtractorOutput.characteristic())
                .actionableEvents(actionableEvents)
                .build();
    }

    @NotNull
    private static ExtractionResult toExtractionResult(@NotNull Set<ActionableEvent> actionableEvents,
            @NotNull CkbExtractorResult ckbExtractorResult) {
        Set<ActionableHotspot> actionableHotspots = Sets.newHashSet();
        Set<ActionableRange> actionableRanges = Sets.newHashSet();
        Set<ActionableGene> actionableGenes = Sets.newHashSet();
        Set<ActionableFusion> actionableFusions = Sets.newHashSet();
        Set<ActionableCharacteristic> actionableCharacteristics = Sets.newHashSet();

        for (ActionableEvent event : actionableEvents) {

            actionableHotspots.addAll(ActionableEventFactory.toActionableHotspots(event, ckbExtractorResult.hotspots()));
            actionableRanges.addAll(ActionableEventFactory.toActionableRanges(event, ckbExtractorResult.codons()));
            actionableRanges.addAll(ActionableEventFactory.toActionableRanges(event, ckbExtractorResult.exons()));

            if (ckbExtractorResult.geneLevelEvent() != null) {
                actionableGenes.add(ActionableEventFactory.geneLevelEventToActionableGene(event, ckbExtractorResult.geneLevelEvent()));
            }

            if (ckbExtractorResult.knownCopyNumber() != null) {
                actionableGenes.add(ActionableEventFactory.copyNumberToActionableGene(event, ckbExtractorResult.knownCopyNumber()));
            }

            if (ckbExtractorResult.knownFusionPair() != null) {
                actionableFusions.add(ActionableEventFactory.toActionableFusion(event, ckbExtractorResult.knownFusionPair()));
            }

            if (ckbExtractorResult.characteristic() != null) {
                actionableCharacteristics.add(ActionableEventFactory.toActionableCharacteristic(event,
                        ckbExtractorResult.characteristic()));
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