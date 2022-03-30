package com.hartwig.hmftools.serve.sources.iclusion;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.iclusion.datamodel.IclusionMutation;
import com.hartwig.hmftools.iclusion.datamodel.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
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
import com.hartwig.hmftools.serve.util.ProgressTracker;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class IclusionExtractor {

    private static final Logger LOGGER = LogManager.getLogger(IclusionExtractor.class);

    @NotNull
    private final EventExtractor eventExtractor;
    @NotNull
    private final ActionableTrialFactory actionableTrialFactory;

    IclusionExtractor(@NotNull final EventExtractor eventExtractor, @NotNull final ActionableTrialFactory actionableTrialFactory) {
        this.eventExtractor = eventExtractor;
        this.actionableTrialFactory = actionableTrialFactory;
    }

    @NotNull
    public ExtractionResult extract(@NotNull List<IclusionTrial> trials) {
        // We assume filtered trials (no empty acronyms, only OR mutations, and no negated mutations
        List<ExtractionResult> extractions = Lists.newArrayList();

        ProgressTracker tracker = new ProgressTracker("iClusion", trials.size());
        for (IclusionTrial trial : trials) {
            for (IclusionMutationCondition mutationCondition : trial.mutationConditions()) {
                for (IclusionMutation mutation : mutationCondition.mutations()) {
                    LOGGER.debug("Interpreting '{}' on '{}' for {}", mutation.name(), mutation.gene(), trial.acronym());

                    if (mutation.type() == EventType.UNKNOWN) {
                        LOGGER.warn("No event type known for '{}' on '{}'", mutation.name(), mutation.gene());
                    } else {
                        extractions.add(extract(trial, mutation));
                    }
                }
            }

            tracker.update();
        }

        return ExtractionFunctions.merge(extractions);
    }

    @NotNull
    private ExtractionResult extract(@NotNull IclusionTrial trial, @NotNull IclusionMutation mutation) {
        String sourceEvent = mutation.gene() + " " + mutation.name();
        List<ActionableTrial> actionableTrials = actionableTrialFactory.toActionableTrials(trial, sourceEvent);
        for (ActionableTrial actionableTrial : actionableTrials) {
            LOGGER.debug("Generated {} based off {}", actionableTrial, trial);
        }

        EventExtractorOutput extraction = eventExtractor.extract(mutation.gene(), null, mutation.type(), mutation.name());

        EventInterpretation interpretation = ImmutableEventInterpretation.builder()
                .source(Knowledgebase.ICLUSION)
                .sourceEvent(sourceEvent)
                .interpretedGene(mutation.gene())
                .interpretedEvent(mutation.name())
                .interpretedEventType(mutation.type())
                .build();

        return toExtractionResult(actionableTrials, extraction, interpretation);
    }

    @NotNull
    private static ExtractionResult toExtractionResult(@NotNull List<ActionableTrial> actionableTrials,
            @NotNull EventExtractorOutput extraction, @NotNull EventInterpretation interpretation) {
        Set<ActionableHotspot> actionableHotspots = Sets.newHashSet();
        Set<ActionableRange> actionableRanges = Sets.newHashSet();
        Set<ActionableGene> actionableGenes = Sets.newHashSet();
        Set<ActionableFusion> actionableFusions = Sets.newHashSet();
        Set<ActionableCharacteristic> actionableCharacteristics = Sets.newHashSet();
        Set<ActionableHLA> actionableHLA = Sets.newHashSet();

        for (ActionableTrial trial : actionableTrials) {
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

            if (extraction.hla() != null) {
                actionableHLA.add(ActionableEventFactory.toActionableHLa(trial, extraction.hla()));
            }
        }

        return ImmutableExtractionResult.builder()
                .refGenomeVersion(Knowledgebase.ICLUSION.refGenomeVersion())
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