package com.hartwig.hmftools.serve.sources.iclusion;

import java.util.HashSet;
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

                    List<ActionableTrial> actionableTrials =
                            actionableTrialFactory.toActionableTrials(trial, mutation.gene() + " " + mutation.name());
                    for (ActionableTrial actionableTrial : actionableTrials) {
                        LOGGER.debug("Generated {} based off {}", actionableTrial, trial);
                    }

                    LOGGER.debug("Interpreting '{}' on '{}' for {}", mutation.name(), mutation.gene(), trial.acronym());
                    if (mutation.type() == EventType.UNKNOWN) {
                        LOGGER.warn("No event type known for '{}' on '{}'", mutation.name(), mutation.gene());
                    }
                    EventExtractorOutput eventExtractorOutput =
                            eventExtractor.extract(mutation.gene(), null, mutation.type(), mutation.name());

                    EventInterpretation interpretation = ImmutableEventInterpretation.builder()
                            .source(Knowledgebase.ICLUSION)
                            .sourceEvent(mutation.name())
                            .interpretedGene(mutation.gene())
                            .interpretedEvent(mutation.name())
                            .interpretedEventType(mutation.type())
                            .build();
                    Set<ActionableTrial> actionableTrialsUniq = new HashSet<>(actionableTrials);
                    extractions.add(toExtractionResult(actionableTrialsUniq, eventExtractorOutput, interpretation));
                }
            }

            tracker.update();
        }

        return ExtractionFunctions.merge(extractions);
    }

    @NotNull
    private static ExtractionResult toExtractionResult(@NotNull Set<ActionableTrial> actionableTrials,
            @NotNull EventExtractorOutput eventExtractions, @NotNull EventInterpretation interpretation) {
        Set<ActionableHotspot> actionableHotspots = Sets.newHashSet();
        Set<ActionableRange> actionableRanges = Sets.newHashSet();
        Set<ActionableGene> actionableGenes = Sets.newHashSet();
        Set<ActionableFusion> actionableFusions = Sets.newHashSet();
        Set<ActionableCharacteristic> actionableCharacteristics = Sets.newHashSet();
        Set<ActionableHLA> actionableHLA = Sets.newHashSet();

        for (ActionableTrial trial : actionableTrials) {

            actionableHotspots.addAll(ActionableEventFactory.toActionableHotspots(trial, eventExtractions.hotspots()));
            actionableRanges.addAll(ActionableEventFactory.toActionableRanges(trial, eventExtractions.codons()));
            actionableRanges.addAll(ActionableEventFactory.toActionableRanges(trial, eventExtractions.exons()));

            if (eventExtractions.geneLevelEvent() != null) {
                actionableGenes.add(ActionableEventFactory.geneLevelEventToActionableGene(trial, eventExtractions.geneLevelEvent()));
            }

            if (eventExtractions.knownCopyNumber() != null) {
                actionableGenes.add(ActionableEventFactory.copyNumberToActionableGene(trial, eventExtractions.knownCopyNumber()));
            }

            if (eventExtractions.knownFusionPair() != null) {
                actionableFusions.add(ActionableEventFactory.toActionableFusion(trial, eventExtractions.knownFusionPair()));
            }

            if (eventExtractions.characteristic() != null) {
                actionableCharacteristics.add(ActionableEventFactory.toActionableCharacteristic(trial, eventExtractions.characteristic()));
            }

            if (eventExtractions.hla() != null) {
                actionableHLA.add(ActionableEventFactory.toActionableHLa(trial, eventExtractions.hla()));
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