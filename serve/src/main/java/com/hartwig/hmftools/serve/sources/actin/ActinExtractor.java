package com.hartwig.hmftools.serve.sources.actin;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.immuno.ActionableHLA;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.extraction.ActionableEventFactory;
import com.hartwig.hmftools.serve.extraction.EventExtractor;
import com.hartwig.hmftools.serve.extraction.EventExtractorOutput;
import com.hartwig.hmftools.serve.extraction.ExtractionFunctions;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ImmutableExtractionResult;
import com.hartwig.hmftools.serve.extraction.events.EventInterpretation;
import com.hartwig.hmftools.serve.extraction.events.ImmutableEventInterpretation;
import com.hartwig.hmftools.serve.sources.actin.classification.ActinEventExtractor;
import com.hartwig.hmftools.serve.sources.actin.classification.ActinEventTypeExtractor;
import com.hartwig.hmftools.serve.sources.actin.classification.ActinKeywords;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.util.ProgressTracker;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ActinExtractor {

    private static final Logger LOGGER = LogManager.getLogger(ActinExtractor.class);

    @NotNull
    private final EventExtractor eventExtractor;

    ActinExtractor(@NotNull final EventExtractor eventExtractor) {
        this.eventExtractor = eventExtractor;
    }

    @NotNull
    public ExtractionResult extract(@NotNull List<ActinEntry> entries) {
        List<ExtractionResult> extractions = Lists.newArrayList();

        ProgressTracker tracker = new ProgressTracker("ACTIN", entries.size());
        for (ActinEntry entry : entries) {
            Set<String> events = ActinEventExtractor.extractEvents(entry);

            for (String event : events) {
                EventType type = ActinEventTypeExtractor.determineEventType(entry, event);
                if (type == EventType.UNKNOWN) {
                    LOGGER.warn("No event type known for '{}'", entry);
                } else {
                    extractions.add(toExtractionResult(entry, event, type));
                }
            }

            tracker.update();
        }

        return ExtractionFunctions.merge(extractions);
    }

    @NotNull
    private ExtractionResult toExtractionResult(@NotNull ActinEntry entry, @NotNull String event, @NotNull EventType type) {
        String gene;
        String sourceEvent;

        if (entry.gene() != null) {
            gene = entry.gene();
            sourceEvent = entry.rule() + ": " + entry.gene();
        } else {
            gene = ActinKeywords.NO_GENE;
            sourceEvent = entry.rule().toString();
        }

        if (entry.mutation() != null) {
            sourceEvent += (" " + entry.mutation());
        }

        EventExtractorOutput extraction = eventExtractor.extract(gene, null, type, event);
        EventInterpretation interpretation = ImmutableEventInterpretation.builder()
                .source(Knowledgebase.ACTIN)
                .sourceEvent(sourceEvent)
                .interpretedGene(gene)
                .interpretedEvent(event)
                .interpretedEventType(type)
                .build();

        return toExtractionResult(ActinTrialFactory.toActinTrial(entry, sourceEvent), extraction, interpretation);
    }

    @NotNull
    private static ExtractionResult toExtractionResult(@NotNull ActinTrial trial, @NotNull EventExtractorOutput extraction,
            @NotNull EventInterpretation interpretation) {
        Set<ActionableHotspot> actionableHotspots = ActionableEventFactory.toActionableHotspots(trial, extraction.hotspots());

        Set<ActionableRange> actionableRanges = Sets.newHashSet();
        actionableRanges.addAll(ActionableEventFactory.toActionableRanges(trial, extraction.codons()));
        actionableRanges.addAll(ActionableEventFactory.toActionableRanges(trial, extraction.exons()));

        Set<ActionableGene> actionableGenes = Sets.newHashSet();
        if (extraction.geneLevelEvent() != null) {
            actionableGenes.add(ActionableEventFactory.geneLevelEventToActionableGene(trial, extraction.geneLevelEvent()));
        }

        if (extraction.knownCopyNumber() != null) {
            actionableGenes.add(ActionableEventFactory.copyNumberToActionableGene(trial, extraction.knownCopyNumber()));
        }

        Set<ActionableFusion> actionableFusions = Sets.newHashSet();
        if (extraction.knownFusionPair() != null) {
            actionableFusions.add(ActionableEventFactory.toActionableFusion(trial, extraction.knownFusionPair()));
        }

        Set<ActionableCharacteristic> actionableCharacteristics = Sets.newHashSet();
        if (extraction.characteristic() != null) {
            actionableCharacteristics.add(ActionableEventFactory.toActionableCharacteristic(trial, extraction.characteristic()));
        }

        Set<ActionableHLA> actionableHLA = Sets.newHashSet();
        if (extraction.hla() != null) {
            actionableHLA.add(ActionableEventFactory.toActionableHLa(trial, extraction.hla()));
        }

        return ImmutableExtractionResult.builder()
                .refGenomeVersion(Knowledgebase.ACTIN.refGenomeVersion())
                .addEventInterpretations((interpretation))
                .actionableHotspots(actionableHotspots)
                .actionableRanges(actionableRanges)
                .actionableGenes(actionableGenes)
                .actionableFusions(actionableFusions)
                .actionableCharacteristics(actionableCharacteristics)
                .actionableHLA(actionableHLA)
                .build();
    }
}