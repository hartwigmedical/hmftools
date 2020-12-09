package com.hartwig.hmftools.serve.sources.vicc;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.fusion.ImmutableActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.hotspot.ImmutableActionableHotspot;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.actionability.range.ImmutableActionableRange;
import com.hartwig.hmftools.serve.actionability.signature.ActionableSignature;
import com.hartwig.hmftools.serve.actionability.signature.ImmutableActionableSignature;
import com.hartwig.hmftools.serve.extraction.EventExtractor;
import com.hartwig.hmftools.serve.extraction.EventExtractorOutput;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ImmutableExtractionResult;
import com.hartwig.hmftools.serve.extraction.codon.CodonAnnotation;
import com.hartwig.hmftools.serve.extraction.copynumber.CopyNumberFunctions;
import com.hartwig.hmftools.serve.extraction.copynumber.ImmutableKnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.exon.ExonAnnotation;
import com.hartwig.hmftools.serve.extraction.fusion.FusionFunctions;
import com.hartwig.hmftools.serve.extraction.fusion.ImmutableKnownFusionPair;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.extraction.hotspot.HotspotFunctions;
import com.hartwig.hmftools.serve.extraction.hotspot.ImmutableKnownHotspot;
import com.hartwig.hmftools.serve.extraction.hotspot.KnownHotspot;
import com.hartwig.hmftools.serve.extraction.signature.SignatureName;
import com.hartwig.hmftools.serve.util.ProgressTracker;
import com.hartwig.hmftools.vicc.annotation.ProteinAnnotationExtractor;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ViccExtractor {

    private static final Logger LOGGER = LogManager.getLogger(ViccExtractor.class);

    @NotNull
    private final EventExtractor eventExtractor;
    @Nullable
    private final String featureInterpretationTsv;

    public ViccExtractor(@NotNull final EventExtractor eventExtractor, @Nullable final String featureInterpretationTsv) {
        this.eventExtractor = eventExtractor;
        this.featureInterpretationTsv = featureInterpretationTsv;
    }

    @NotNull
    public ExtractionResult extract(@NotNull List<ViccEntry> entries) throws IOException {
        Map<ViccEntry, ViccExtractionResult> resultsPerEntry = Maps.newHashMap();

        ProgressTracker tracker = new ProgressTracker("VICC", entries.size());
        for (ViccEntry entry : entries) {
            resultsPerEntry.put(entry, extractSingleEntry(entry));

            tracker.update();
        }

        ViccUtil.printExtractionResults(resultsPerEntry);

        if (featureInterpretationTsv != null) {
            ViccUtil.writeInterpretationToTsv(featureInterpretationTsv, resultsPerEntry);
        }

        ImmutableExtractionResult.Builder outputBuilder = ImmutableExtractionResult.builder()
                .knownHotspots(convertToHotspots(resultsPerEntry))
                .knownCopyNumbers(convertToKnownAmpsDels(resultsPerEntry))
                .knownFusionPairs(convertToKnownFusions(resultsPerEntry));

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
        Map<Feature, SignatureName> signaturesPerFeature = Maps.newHashMap();

        for (Feature feature : entry.features()) {
            String gene = feature.geneSymbol();
            if (gene == null) {
                LOGGER.warn("No gene configured for {}. Skipping!", feature);
            } else {
                EventExtractorOutput extractorOutput = eventExtractor.extract(gene, entry.transcriptId(), feature.type(), feature.name());
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

                if (extractorOutput.signatureName() != null) {
                    signaturesPerFeature.put(feature, extractorOutput.signatureName());
                }
            }
        }

        ActionableEvent actionableEvidence = ActionableEvidenceFactory.toActionableEvent(entry);

        return ImmutableViccExtractionResult.builder()
                .hotspotsPerFeature(hotspotsPerFeature)
                .codonsPerFeature(codonsPerFeature)
                .exonsPerFeature(exonsPerFeature)
                .geneLevelEventsPerFeature(geneLevelEventsPerFeature)
                .ampsDelsPerFeature(ampsDelsPerFeature)
                .fusionsPerFeature(fusionsPerFeature)
                .signaturesPerFeature(signaturesPerFeature)
                .actionableEvent(actionableEvidence)
                .build();
    }

    @NotNull
    private static Set<KnownHotspot> convertToHotspots(@NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
        ProteinAnnotationExtractor proteinExtractor = new ProteinAnnotationExtractor();
        Set<KnownHotspot> hotspots = Sets.newHashSet();
        for (Map.Entry<ViccEntry, ViccExtractionResult> entryResult : resultsPerEntry.entrySet()) {
            ViccEntry entry = entryResult.getKey();
            Knowledgebase source = ViccSource.toKnowledgebase(entry.source());
            for (Map.Entry<Feature, List<VariantHotspot>> featureResult : entryResult.getValue().hotspotsPerFeature().entrySet()) {
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
        }

        return HotspotFunctions.consolidate(hotspots);
    }

    @NotNull
    private static Set<KnownCopyNumber> convertToKnownAmpsDels(@NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
        Set<KnownCopyNumber> copyNumbers = Sets.newHashSet();
        for (Map.Entry<ViccEntry, ViccExtractionResult> entry : resultsPerEntry.entrySet()) {
            Knowledgebase source = ViccSource.toKnowledgebase(entry.getKey().source());
            for (KnownCopyNumber copyNumber : entry.getValue().ampsDelsPerFeature().values()) {
                copyNumbers.add(ImmutableKnownCopyNumber.builder().from(copyNumber).addSources(source).build());
            }
        }

        return CopyNumberFunctions.consolidate(copyNumbers);
    }

    @NotNull
    private static Set<KnownFusionPair> convertToKnownFusions(@NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
        Set<KnownFusionPair> fusions = Sets.newHashSet();
        for (Map.Entry<ViccEntry, ViccExtractionResult> entry : resultsPerEntry.entrySet()) {
            Knowledgebase source = ViccSource.toKnowledgebase(entry.getKey().source());
            for (KnownFusionPair fusionPair : entry.getValue().fusionsPerFeature().values()) {
                fusions.add(ImmutableKnownFusionPair.builder().from(fusionPair).addSources(source).build());
            }
        }

        return FusionFunctions.consolidate(fusions);
    }

    private static void addActionability(@NotNull ImmutableExtractionResult.Builder outputBuilder,
            @NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
        Set<ActionableHotspot> actionableHotspots = Sets.newHashSet();
        Set<ActionableRange> actionableRanges = Sets.newHashSet();
        Set<ActionableGene> actionableGenes = Sets.newHashSet();
        Set<ActionableFusion> actionableFusions = Sets.newHashSet();
        Set<ActionableSignature> actionableSignatures = Sets.newHashSet();

        for (Map.Entry<ViccEntry, ViccExtractionResult> entry : resultsPerEntry.entrySet()) {
            ViccExtractionResult result = entry.getValue();
            ActionableEvent event = result.actionableEvent();
            if (event != null) {
                actionableHotspots.addAll(extractActionableHotspots(event, result.hotspotsPerFeature().values()));
                actionableRanges.addAll(extractActionableRanges(event,
                        result.codonsPerFeature().values(),
                        result.exonsPerFeature().values()));
                actionableGenes.addAll(extractActionableAmpsDels(event, result.ampsDelsPerFeature().values()));
                actionableGenes.addAll(extractActionableGeneLevelEvents(event, result.geneLevelEventsPerFeature().values()));
                actionableFusions.addAll(extractActionableFusions(event, result.fusionsPerFeature().values()));
                actionableSignatures.addAll(extractActionableSignatures(event, result.signaturesPerFeature().values()));
            }
        }

        outputBuilder.actionableHotspots(actionableHotspots);
        outputBuilder.actionableRanges(actionableRanges);
        outputBuilder.actionableGenes(actionableGenes);
        outputBuilder.actionableFusions(actionableFusions);
        outputBuilder.actionableSignatures(actionableSignatures);
    }

    @NotNull
    private static Set<ActionableHotspot> extractActionableHotspots(@NotNull ActionableEvent actionableEvent,
            @NotNull Iterable<List<VariantHotspot>> hotspotLists) {
        Set<ActionableHotspot> actionableHotspots = Sets.newHashSet();
        for (List<VariantHotspot> hotspotList : hotspotLists) {
            for (VariantHotspot hotspot : hotspotList) {
                actionableHotspots.add(ImmutableActionableHotspot.builder()
                        .from(actionableEvent)
                        .chromosome(hotspot.chromosome())
                        .position(hotspot.position())
                        .ref(hotspot.ref())
                        .alt(hotspot.alt())
                        .build());
            }
        }
        return actionableHotspots;
    }

    @NotNull
    private static Set<ActionableRange> extractActionableRanges(@NotNull ActionableEvent actionableEvent,
            @NotNull Iterable<List<CodonAnnotation>> codonLists, @NotNull Iterable<List<ExonAnnotation>> exonLists) {
        Set<ActionableRange> actionableRanges = Sets.newHashSet();
        for (List<CodonAnnotation> codonList : codonLists) {
            for (CodonAnnotation codon : codonList) {
                actionableRanges.add(ImmutableActionableRange.builder()
                        .from(actionableEvent)
                        .chromosome(codon.chromosome())
                        .start(codon.start())
                        .end(codon.end())
                        .gene(codon.gene())
                        .mutationType(codon.mutationType())
                        .build());
            }
        }

        for (List<ExonAnnotation> exonList : exonLists) {
            for (ExonAnnotation exon : exonList) {
                actionableRanges.add(ImmutableActionableRange.builder()
                        .from(actionableEvent)
                        .chromosome(exon.chromosome())
                        .start(exon.start())
                        .end(exon.end())
                        .gene(exon.gene())
                        .mutationType(exon.mutationType())
                        .build());
            }
        }
        return actionableRanges;
    }

    @NotNull
    private static Set<ActionableGene> extractActionableAmpsDels(@NotNull ActionableEvent actionableEvent,
            @NotNull Iterable<KnownCopyNumber> copyNumbers) {
        Set<ActionableGene> actionableGenes = Sets.newHashSet();
        for (KnownCopyNumber copyNumber : copyNumbers) {
            GeneLevelEvent event;
            switch (copyNumber.type()) {
                case AMPLIFICATION: {
                    event = GeneLevelEvent.AMPLIFICATION;
                    break;
                }
                case DELETION: {
                    event = GeneLevelEvent.DELETION;
                    break;
                }
                default:
                    throw new IllegalStateException("Invalid copy number type: " + copyNumber.type());
            }

            actionableGenes.add(ImmutableActionableGene.builder().from(actionableEvent).gene(copyNumber.gene()).event(event).build());
        }
        return actionableGenes;
    }

    @NotNull
    private static Set<ActionableGene> extractActionableGeneLevelEvents(@NotNull ActionableEvent actionableEvent,
            @NotNull Iterable<GeneLevelAnnotation> geneLevelEvents) {
        Set<ActionableGene> actionableGenes = Sets.newHashSet();
        for (GeneLevelAnnotation geneLevelEvent : geneLevelEvents) {
            actionableGenes.add(ImmutableActionableGene.builder()
                    .from(actionableEvent)
                    .gene(geneLevelEvent.gene())
                    .event(geneLevelEvent.event())
                    .build());
        }
        return actionableGenes;
    }

    @NotNull
    private static Set<ActionableFusion> extractActionableFusions(@NotNull ActionableEvent actionableEvent,
            @NotNull Iterable<KnownFusionPair> fusionAnnotations) {
        Set<ActionableFusion> actionableFusions = Sets.newHashSet();
        for (KnownFusionPair fusion : fusionAnnotations) {
            actionableFusions.add(ImmutableActionableFusion.builder()
                    .from(actionableEvent)
                    .geneUp(fusion.geneUp())
                    .minExonUp(fusion.minExonUp())
                    .maxExonUp(fusion.maxExonUp())
                    .geneDown(fusion.geneDown())
                    .minExonDown(fusion.minExonDown())
                    .maxExonDown(fusion.maxExonDown())
                    .build());
        }
        return actionableFusions;
    }

    @NotNull
    private static Set<ActionableSignature> extractActionableSignatures(@NotNull ActionableEvent actionableEvent,
            @NotNull Iterable<SignatureName> signatures) {
        Set<ActionableSignature> actionableSignatures = Sets.newHashSet();
        for (SignatureName signature : signatures) {
            actionableSignatures.add(ImmutableActionableSignature.builder().from(actionableEvent).name(signature).build());
        }
        return actionableSignatures;
    }
}
