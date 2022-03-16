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
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.ckb.CkbExtractor;
import com.hartwig.hmftools.serve.util.ProgressTracker;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ActinExtractor {

    private static final Logger LOGGER = LogManager.getLogger(CkbExtractor.class);

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
            ActinTrial trial = ActinTrialFactory.toActinTrial(entry, entry.rule().name());

            for (String event : events) {
                String eventType = event;
                if (event.split(" ").length > 1) {
                    eventType = event.split(" ", 2)[0];
                }

                EventType type = ActinEventTypeExtractor.determineEventType(entry, eventType);
                if (type == EventType.UNKNOWN) {
                    LOGGER.warn("No event type known for '{}' on '{}'", entry, entry.gene());
                } else {
                    String gene = entry.gene() != null ? entry.gene() : "-";
                    EventInterpretation interpretation = ImmutableEventInterpretation.builder()
                            .source(Knowledgebase.ACTIN)
                            .sourceEvent(entry.rule().toString())
                            .interpretedGene(gene)
                            .interpretedEvent(event)
                            .interpretedEventType(type)
                            .build();

                    extractions.add(toExtractionResult(trial, eventExtractor.extract(gene, null, type, event), interpretation));
                }
            }

            tracker.update();
        }

        return ExtractionFunctions.merge(extractions);
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
                .eventInterpretations(Lists.newArrayList(interpretation))
                .actionableHotspots(actionableHotspots)
                .actionableRanges(actionableRanges)
                .actionableGenes(actionableGenes)
                .actionableFusions(actionableFusions)
                .actionableCharacteristics(actionableCharacteristics)
                .actionableHLA(actionableHLA)
                .build();
    }
}