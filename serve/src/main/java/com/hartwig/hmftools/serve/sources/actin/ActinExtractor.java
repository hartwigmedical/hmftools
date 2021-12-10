package com.hartwig.hmftools.serve.sources.actin;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.EventClassifierFactory;
import com.hartwig.hmftools.common.serve.classification.EventType;
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
import com.hartwig.hmftools.serve.sources.actin.classification.ActinClassificationConfig;
import com.hartwig.hmftools.serve.sources.actin.classification.ActinEventTypeExtractor;
import com.hartwig.hmftools.serve.util.ProgressTracker;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ActinExtractor {

    private static final Logger LOGGER = LogManager.getLogger(ActinExtractor.class);
    private static final EventClassifier CLASSIFIER = EventClassifierFactory.buildClassifier(ActinClassificationConfig.build());
    private static final ActinEventTypeExtractor EXTRACTOR = new ActinEventTypeExtractor();

    @NotNull
    private final EventExtractor eventExtractor;
    @NotNull
    private final com.hartwig.hmftools.serve.sources.actin.ActionableTrialFactory actionableTrialFactory;

    ActinExtractor(@NotNull final EventExtractor eventExtractor,
            @NotNull final com.hartwig.hmftools.serve.sources.actin.ActionableTrialFactory actionableTrialFactory) {
        this.eventExtractor = eventExtractor;
        this.actionableTrialFactory = actionableTrialFactory;
    }

    @NotNull
    public ExtractionResult extract(@NotNull List<ActinTrial> trials) {
        // We assume filtered trials (no empty acronyms, only OR mutations, and no negated mutations

        ProgressTracker tracker = new ProgressTracker("ACTIN", trials.size());
        List<ExtractionResult> extractions = Lists.newArrayList();
        for (ActinTrial trial : trials) {
            List<com.hartwig.hmftools.serve.sources.actin.ActionableTrial> actionableTrials =
                    actionableTrialFactory.toActionableEntries(trial);

            List<EventExtractorOutput> eventExtractions = Lists.newArrayList();

            String gene = EXTRACTOR.extractGene(trials);
            String event = EXTRACTOR.extractEvent(trials);

            EventType eventType = CLASSIFIER.determineType(gene, event);

            eventExtractions.add(eventExtractor.extract(gene, null, eventType, event));

            extractions.add(toExtractionResult(actionableTrials, eventExtractions));

            tracker.update();
        }

        return ExtractionFunctions.merge(extractions);
    }

    @NotNull
    private static ExtractionResult toExtractionResult(
            @NotNull List<com.hartwig.hmftools.serve.sources.actin.ActionableTrial> actionableTrials,
            @NotNull List<EventExtractorOutput> eventExtractions) {
        Set<ActionableHotspot> actionableHotspots = Sets.newHashSet();
        Set<ActionableRange> actionableRanges = Sets.newHashSet();
        Set<ActionableGene> actionableGenes = Sets.newHashSet();
        Set<ActionableFusion> actionableFusions = Sets.newHashSet();
        Set<ActionableCharacteristic> actionableCharacteristics = Sets.newHashSet();

        for (ActionableTrial trial : actionableTrials) {
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
                .refGenomeVersion(Knowledgebase.ACTIN.refGenomeVersion())
                .actionableHotspots(actionableHotspots)
                .actionableRanges(actionableRanges)
                .actionableGenes(actionableGenes)
                .actionableFusions(actionableFusions)
                .actionableCharacteristics(actionableCharacteristics)
                .build();
    }
}
