package com.hartwig.hmftools.serve.sources.vicc;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
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
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristic;
import com.hartwig.hmftools.serve.extraction.codon.CodonAnnotation;
import com.hartwig.hmftools.serve.extraction.codon.CodonFunctions;
import com.hartwig.hmftools.serve.extraction.codon.ImmutableKnownCodon;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodon;
import com.hartwig.hmftools.serve.extraction.copynumber.CopyNumberFunctions;
import com.hartwig.hmftools.serve.extraction.copynumber.ImmutableKnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.events.EventInterpretation;
import com.hartwig.hmftools.serve.extraction.events.ImmutableEventInterpretation;
import com.hartwig.hmftools.serve.extraction.exon.ExonAnnotation;
import com.hartwig.hmftools.serve.extraction.exon.ExonFunctions;
import com.hartwig.hmftools.serve.extraction.exon.ImmutableKnownExon;
import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.fusion.FusionFunctions;
import com.hartwig.hmftools.serve.extraction.fusion.ImmutableKnownFusionPair;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.extraction.hotspot.HotspotFunctions;
import com.hartwig.hmftools.serve.extraction.hotspot.ImmutableKnownHotspot;
import com.hartwig.hmftools.serve.extraction.hotspot.KnownHotspot;
import com.hartwig.hmftools.serve.extraction.immuno.ImmunoHLA;
import com.hartwig.hmftools.serve.extraction.range.RangeAnnotation;
import com.hartwig.hmftools.serve.util.ProgressTracker;
import com.hartwig.hmftools.vicc.annotation.ViccProteinAnnotationExtractor;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ViccExtractor {

    private static final Logger LOGGER = LogManager.getLogger(ViccExtractor.class);

    @NotNull
    private final EventExtractor eventExtractor;
    @NotNull
    private final ActionableEvidenceFactory actionableEvidenceFactory;

    public ViccExtractor(@NotNull final EventExtractor eventExtractor, @NotNull final ActionableEvidenceFactory actionableEvidenceFactory) {
        this.eventExtractor = eventExtractor;
        this.actionableEvidenceFactory = actionableEvidenceFactory;
    }

    @NotNull
    public ExtractionResult extract(@NotNull List<ViccEntry> entries) {
        Map<ViccEntry, ViccExtractionResult> resultsPerEntry = Maps.newHashMap();
        List<ExtractionResult> extractions = Lists.newArrayList();

        ProgressTracker tracker = new ProgressTracker("VICC", entries.size());
        for (ViccEntry entry : entries) {
            resultsPerEntry.put(entry, extractSingleEntry(entry));
            tracker.update();
        }

        for (Map.Entry<ViccEntry, ViccExtractionResult> entryResult : resultsPerEntry.entrySet()) {
            extractions.add(toExtractionResult(entryResult.getValue(), entryResult.getKey()));
        }

        actionableEvidenceFactory.evaluateCuration();
        ViccUtil.printExtractionResults(resultsPerEntry);
        return ExtractionFunctions.merge(extractions);
    }

    @NotNull
    private static ExtractionResult toExtractionResult(ViccExtractionResult resultsPerEntry, @NotNull ViccEntry entry) {
        // Assume all VICC knowledgebases are on the same ref genome version
        ImmutableExtractionResult.Builder outputBuilder = ImmutableExtractionResult.builder()
                .eventInterpretations(Lists.newArrayList(resultsPerEntry.eventInterpretations()))
                .refGenomeVersion(Knowledgebase.VICC_CGI.refGenomeVersion())
                .knownHotspots(convertToHotspots(resultsPerEntry, entry))
                .knownCodons(convertToCodons(resultsPerEntry, entry))
                .knownExons(convertToExons(resultsPerEntry, entry))
                .knownCopyNumbers(convertToKnownAmpsDels(resultsPerEntry, entry))
                .knownFusionPairs(convertToKnownFusions(resultsPerEntry, entry));

        addActionability(outputBuilder, resultsPerEntry);

        return outputBuilder.build();
    }

    @NotNull
    private ViccExtractionResult extractSingleEntry(@NotNull ViccEntry entry) {
        Map<Feature, List<VariantHotspot>> hotspotsPerFeature = Maps.newHashMap();
        Map<Feature, List<CodonAnnotation>> codonsPerFeature = Maps.newHashMap();
        Map<Feature, List<ExonAnnotation>> exonsPerFeature = Maps.newHashMap();
        Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = Maps.newHashMap();
        Map<Feature, KnownCopyNumber> ampsDelsPerFeature = Maps.newHashMap();
        Map<Feature, KnownFusionPair> fusionsPerFeature = Maps.newHashMap();
        Map<Feature, TumorCharacteristic> characteristicsPerFeature = Maps.newHashMap();
        Map<Feature, ImmunoHLA> hlaPerFeature = Maps.newHashMap();
        List<EventInterpretation> eventInterpretations = Lists.newArrayList();
        ImmutableViccExtractionResult.Builder viccBuilderBuilder = ImmutableViccExtractionResult.builder();

        for (Feature feature : entry.features()) {
            String gene = feature.geneSymbol();
            if (gene == null) {
                LOGGER.warn("No gene configured for {}. Skipping!", feature);
            } else {
                EventExtractorOutput extractorOutput = eventExtractor.extract(gene, entry.transcriptId(), feature.type(), feature.name());

                eventInterpretations.add(ImmutableEventInterpretation.builder()
                        .source(ActionableEvidenceFactory.fromViccSource(entry.source()))
                        .sourceEvent(feature.name())
                        .interpretedGene(gene)
                        .interpretedEvent(feature.name())
                        .interpretedEventType(feature.type())
                        .build());

                //TODO: sourceEvent will be empty in v1.9. Will be fixed later.
                String sourceEvent = gene + " " + feature.name();
                Set<ActionableEvent> actionableEvents = actionableEvidenceFactory.toActionableEvents(entry, sourceEvent);

                if (extractorOutput.hotspots() != null) {
                    hotspotsPerFeature.put(feature, extractorOutput.hotspots());
                }

                if (extractorOutput.codons() != null) {
                    codonsPerFeature.put(feature, extractorOutput.codons());
                }

                if (extractorOutput.exons() != null) {
                    exonsPerFeature.put(feature, extractorOutput.exons());
                }

                if (extractorOutput.geneLevelEvent() != null) {
                    geneLevelEventsPerFeature.put(feature, extractorOutput.geneLevelEvent());
                }

                if (extractorOutput.knownCopyNumber() != null) {
                    ampsDelsPerFeature.put(feature, extractorOutput.knownCopyNumber());
                }

                if (extractorOutput.knownFusionPair() != null) {
                    fusionsPerFeature.put(feature, extractorOutput.knownFusionPair());
                }

                if (extractorOutput.characteristic() != null) {
                    characteristicsPerFeature.put(feature, extractorOutput.characteristic());
                }

                if (extractorOutput.hla() != null) {
                    hlaPerFeature.put(feature, extractorOutput.hla());
                }

                viccBuilderBuilder.actionableEvents(actionableEvents);

            }
        }
        viccBuilderBuilder.eventInterpretations(eventInterpretations);
        viccBuilderBuilder.hotspotsPerFeature(hotspotsPerFeature);
        viccBuilderBuilder.codonsPerFeature(codonsPerFeature);
        viccBuilderBuilder.exonsPerFeature(exonsPerFeature);
        viccBuilderBuilder.geneLevelEventsPerFeature(geneLevelEventsPerFeature);
        viccBuilderBuilder.ampsDelsPerFeature(ampsDelsPerFeature);
        viccBuilderBuilder.fusionsPerFeature(fusionsPerFeature);
        viccBuilderBuilder.characteristicsPerFeature(characteristicsPerFeature);
        viccBuilderBuilder.HLAPerFeature(hlaPerFeature);

        return viccBuilderBuilder.build();
    }

    @NotNull
    private static Set<KnownHotspot> convertToHotspots(@NotNull ViccExtractionResult resultsPerEntry, @NotNull ViccEntry entry) {
        ViccProteinAnnotationExtractor proteinExtractor = new ViccProteinAnnotationExtractor();
        Set<KnownHotspot> hotspots = Sets.newHashSet();
        Knowledgebase source = ViccSource.toKnowledgebase(entry.source());
        for (Map.Entry<Feature, List<VariantHotspot>> featureResult : resultsPerEntry.hotspotsPerFeature().entrySet()) {
            Feature feature = featureResult.getKey();
            for (VariantHotspot hotspot : featureResult.getValue()) {
                hotspots.add(ImmutableKnownHotspot.builder()
                        .from(hotspot)
                        .addSources(source)
                        .gene(feature.geneSymbol())
                        .transcript(entry.transcriptId())
                        .proteinAnnotation(proteinExtractor.apply(feature.name()))
                        .build());
            }
        }

        return HotspotFunctions.consolidate(hotspots);
    }

    @NotNull
    private static Set<KnownCodon> convertToCodons(@NotNull ViccExtractionResult resultsPerEntry, @NotNull ViccEntry entry) {
        Set<KnownCodon> codons = Sets.newHashSet();

        Knowledgebase source = ViccSource.toKnowledgebase(entry.source());
        for (List<CodonAnnotation> annotations : resultsPerEntry.codonsPerFeature().values()) {
            for (CodonAnnotation annotation : annotations) {
                codons.add(ImmutableKnownCodon.builder().annotation(annotation).addSources(source).build());
            }
        }

        return CodonFunctions.consolidate(codons);
    }

    @NotNull
    private static Set<KnownExon> convertToExons(@NotNull ViccExtractionResult resultsPerEntry, @NotNull ViccEntry entry) {
        Set<KnownExon> exons = Sets.newHashSet();

        Knowledgebase source = ViccSource.toKnowledgebase(entry.source());
        for (List<ExonAnnotation> annotations : resultsPerEntry.exonsPerFeature().values()) {
            for (ExonAnnotation annotation : annotations) {
                exons.add(ImmutableKnownExon.builder().annotation(annotation).addSources(source).build());
            }
        }

        return ExonFunctions.consolidate(exons);
    }

    @NotNull
    private static Set<KnownCopyNumber> convertToKnownAmpsDels(@NotNull ViccExtractionResult resultsPerEntry, @NotNull ViccEntry entry) {
        Set<KnownCopyNumber> copyNumbers = Sets.newHashSet();
        Knowledgebase source = ViccSource.toKnowledgebase(entry.source());
        for (KnownCopyNumber copyNumber : resultsPerEntry.ampsDelsPerFeature().values()) {
            copyNumbers.add(ImmutableKnownCopyNumber.builder().from(copyNumber).addSources(source).build());
        }

        return CopyNumberFunctions.consolidate(copyNumbers);
    }

    @NotNull
    private static Set<KnownFusionPair> convertToKnownFusions(@NotNull ViccExtractionResult resultsPerEntry, @NotNull ViccEntry entry) {
        Set<KnownFusionPair> fusions = Sets.newHashSet();
        Knowledgebase source = ViccSource.toKnowledgebase(entry.source());
        for (KnownFusionPair fusionPair : resultsPerEntry.fusionsPerFeature().values()) {
            fusions.add(ImmutableKnownFusionPair.builder().from(fusionPair).addSources(source).build());
        }

        return FusionFunctions.consolidate(fusions);
    }

    private static void addActionability(@NotNull ImmutableExtractionResult.Builder outputBuilder, @NotNull ViccExtractionResult result) {
        Set<ActionableHotspot> actionableHotspots = Sets.newHashSet();
        Set<ActionableRange> actionableRanges = Sets.newHashSet();
        Set<ActionableGene> actionableGenes = Sets.newHashSet();
        Set<ActionableFusion> actionableFusions = Sets.newHashSet();
        Set<ActionableCharacteristic> actionableCharacteristics = Sets.newHashSet();

        actionableHotspots.addAll(extractActionableHotspots(result.actionableEvents(), result.hotspotsPerFeature().values()));
        actionableRanges.addAll(extractActionableRanges(result.actionableEvents(), result.codonsPerFeature().values()));
        actionableRanges.addAll(extractActionableRanges(result.actionableEvents(), result.exonsPerFeature().values()));
        actionableGenes.addAll(extractActionableAmpsDels(result.actionableEvents(), result.ampsDelsPerFeature().values()));
        actionableGenes.addAll(extractActionableGeneLevelEvents(result.actionableEvents(), result.geneLevelEventsPerFeature().values()));
        actionableFusions.addAll(extractActionableFusions(result.actionableEvents(), result.fusionsPerFeature().values()));
        actionableCharacteristics.addAll(extractActionableCharacteristics(result.actionableEvents(),
                result.characteristicsPerFeature().values()));

        outputBuilder.actionableHotspots(actionableHotspots);
        outputBuilder.actionableRanges(actionableRanges);
        outputBuilder.actionableGenes(actionableGenes);
        outputBuilder.actionableFusions(actionableFusions);
        outputBuilder.actionableCharacteristics(actionableCharacteristics);
    }

    @NotNull
    private static Set<ActionableHotspot> extractActionableHotspots(@NotNull Iterable<ActionableEvent> actionableEvents,
            @NotNull Iterable<List<VariantHotspot>> hotspotLists) {
        Set<ActionableHotspot> actionableHotspots = Sets.newHashSet();
        for (List<VariantHotspot> hotspotList : hotspotLists) {
            for (ActionableEvent actionableEvent : actionableEvents) {
                actionableHotspots.addAll(ActionableEventFactory.toActionableHotspots(actionableEvent, hotspotList));
            }
        }
        return actionableHotspots;
    }

    @NotNull
    private static <T extends RangeAnnotation> Set<ActionableRange> extractActionableRanges(
            @NotNull Iterable<ActionableEvent> actionableEvents, @NotNull Iterable<List<T>> rangeLists) {
        Set<ActionableRange> actionableRanges = Sets.newHashSet();
        for (List<T> rangeList : rangeLists) {
            for (ActionableEvent actionableEvent : actionableEvents) {
                actionableRanges.addAll(ActionableEventFactory.toActionableRanges(actionableEvent, rangeList));
            }
        }
        return actionableRanges;
    }

    @NotNull
    private static Set<ActionableGene> extractActionableAmpsDels(@NotNull Iterable<ActionableEvent> actionableEvents,
            @NotNull Iterable<KnownCopyNumber> copyNumbers) {
        Set<ActionableGene> actionableGenes = Sets.newHashSet();
        for (KnownCopyNumber copyNumber : copyNumbers) {
            for (ActionableEvent actionableEvent : actionableEvents) {
                actionableGenes.add(ActionableEventFactory.copyNumberToActionableGene(actionableEvent, copyNumber));
            }
        }
        return actionableGenes;
    }

    @NotNull
    private static Set<ActionableGene> extractActionableGeneLevelEvents(@NotNull Iterable<ActionableEvent> actionableEvents,
            @NotNull Iterable<GeneLevelAnnotation> geneLevelEvents) {
        Set<ActionableGene> actionableGenes = Sets.newHashSet();
        for (GeneLevelAnnotation geneLevelEvent : geneLevelEvents) {
            for (ActionableEvent actionableEvent : actionableEvents) {
                actionableGenes.add(ActionableEventFactory.geneLevelEventToActionableGene(actionableEvent, geneLevelEvent));
            }
        }
        return actionableGenes;
    }

    @NotNull
    private static Set<ActionableFusion> extractActionableFusions(@NotNull Iterable<ActionableEvent> actionableEvents,
            @NotNull Iterable<KnownFusionPair> knownFusions) {
        Set<ActionableFusion> actionableFusions = Sets.newHashSet();
        for (KnownFusionPair fusion : knownFusions) {
            for (ActionableEvent actionableEvent : actionableEvents) {
                actionableFusions.add(ActionableEventFactory.toActionableFusion(actionableEvent, fusion));
            }
        }
        return actionableFusions;
    }

    @NotNull
    private static Set<ActionableCharacteristic> extractActionableCharacteristics(@NotNull Iterable<ActionableEvent> actionableEvents,
            @NotNull Iterable<TumorCharacteristic> characteristics) {
        Set<ActionableCharacteristic> actionableCharacteristics = Sets.newHashSet();
        for (TumorCharacteristic characteristic : characteristics) {
            for (ActionableEvent actionableEvent : actionableEvents) {
                actionableCharacteristics.add(ActionableEventFactory.toActionableCharacteristic(actionableEvent, characteristic));
            }
        }
        return actionableCharacteristics;
    }
}