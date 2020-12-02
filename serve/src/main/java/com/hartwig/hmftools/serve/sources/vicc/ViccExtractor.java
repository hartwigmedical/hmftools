package com.hartwig.hmftools.serve.sources.vicc;

import static com.hartwig.hmftools.serve.sources.vicc.ViccUtil.toKnowledgebase;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.serve.classification.MutationType;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.fusion.ImmutableActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.hotspot.ImmutableActionableHotspot;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.actionability.range.ImmutableActionableRange;
import com.hartwig.hmftools.serve.actionability.signature.ActionableSignature;
import com.hartwig.hmftools.serve.actionability.signature.ImmutableActionableSignature;
import com.hartwig.hmftools.serve.actionability.signature.SignatureName;
import com.hartwig.hmftools.serve.copynumber.CopyNumberFunctions;
import com.hartwig.hmftools.serve.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.fusion.FusionFunctions;
import com.hartwig.hmftools.serve.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.hotspot.HotspotFunctions;
import com.hartwig.hmftools.serve.hotspot.ImmutableKnownHotspot;
import com.hartwig.hmftools.serve.hotspot.KnownHotspot;
import com.hartwig.hmftools.serve.sources.ExtractionOutput;
import com.hartwig.hmftools.serve.sources.ImmutableExtractionOutput;
import com.hartwig.hmftools.serve.sources.vicc.annotation.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.annotation.GeneRangeAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.extractor.CopyNumberExtractor;
import com.hartwig.hmftools.serve.sources.vicc.extractor.FusionExtractor;
import com.hartwig.hmftools.serve.sources.vicc.extractor.GeneLevelExtractor;
import com.hartwig.hmftools.serve.sources.vicc.extractor.GeneRangeExtractor;
import com.hartwig.hmftools.serve.sources.vicc.extractor.HotspotExtractor;
import com.hartwig.hmftools.serve.sources.vicc.extractor.SignatureExtractor;
import com.hartwig.hmftools.vicc.annotation.ProteinAnnotationExtractor;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ViccExtractor {

    private static final Logger LOGGER = LogManager.getLogger(ViccExtractor.class);

    @NotNull
    private final HotspotExtractor hotspotExtractor;
    @NotNull
    private final CopyNumberExtractor copyNumberExtractor;
    @NotNull
    private final FusionExtractor fusionExtractor;
    @NotNull
    private final GeneLevelExtractor geneLevelExtractor;
    @NotNull
    private final GeneRangeExtractor geneRangeExtractor;
    @NotNull
    private final SignatureExtractor signatureExtractor;
    @Nullable
    private final String featureInterpretationTsv;

    public ViccExtractor(@NotNull final HotspotExtractor hotspotExtractor, @NotNull final CopyNumberExtractor copyNumberExtractor,
            @NotNull final FusionExtractor fusionExtractor, @NotNull final GeneLevelExtractor geneLevelExtractor,
            @NotNull final GeneRangeExtractor geneRangeExtractor, @NotNull final SignatureExtractor signatureExtractor,
            @Nullable final String featureInterpretationTsv) {
        this.hotspotExtractor = hotspotExtractor;
        this.copyNumberExtractor = copyNumberExtractor;
        this.fusionExtractor = fusionExtractor;
        this.geneLevelExtractor = geneLevelExtractor;
        this.geneRangeExtractor = geneRangeExtractor;
        this.signatureExtractor = signatureExtractor;
        this.featureInterpretationTsv = featureInterpretationTsv;
    }

    @NotNull
    public ExtractionOutput extractFromViccEntries(@NotNull List<ViccEntry> viccEntries) throws IOException {
        Map<ViccEntry, ViccExtractionResult> resultsPerEntry = Maps.newHashMap();

        for (ViccEntry entry : viccEntries) {
            Map<Feature, List<VariantHotspot>> hotspotsPerFeature = hotspotExtractor.extractHotspots(entry);
            Map<Feature, KnownCopyNumber> ampsDelsPerFeature = copyNumberExtractor.extractAmplificationsDeletions(entry);
            Map<Feature, KnownFusionPair> fusionsPerFeature = fusionExtractor.extractFusionPairs(entry);
            Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = geneLevelExtractor.extractGeneLevelEvents(entry);
            Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature = geneRangeExtractor.extractGeneRanges(entry);
            Map<Feature, SignatureName> signaturesPerFeature = signatureExtractor.extractSignatures(entry);

            ActionableEvent actionableEvidence = ActionableEvidenceFactory.toActionableEvent(entry);

            resultsPerEntry.put(entry,
                    ImmutableViccExtractionResult.builder()
                            .hotspotsPerFeature(hotspotsPerFeature)
                            .ampsDelsPerFeature(ampsDelsPerFeature)
                            .fusionsPerFeature(fusionsPerFeature)
                            .geneLevelEventsPerFeature(geneLevelEventsPerFeature)
                            .geneRangesPerFeature(geneRangesPerFeature)
                            .signaturesPerFeature(signaturesPerFeature)
                            .actionableEvent(actionableEvidence)
                            .build());
        }

        printResults(resultsPerEntry);

        if (featureInterpretationTsv != null) {
            writeInterpretationToTsv(resultsPerEntry, featureInterpretationTsv);
        }

        ImmutableExtractionOutput.Builder outputBuilder = ImmutableExtractionOutput.builder()
                .knownHotspots(convertToHotspots(resultsPerEntry, hotspotExtractor.proteinAnnotationExtractor()))
                .knownCopyNumbers(convertToKnownAmpsDels(resultsPerEntry))
                .knownFusionPairs(convertToKnownFusions(resultsPerEntry));

        addActionability(outputBuilder, resultsPerEntry);

        return outputBuilder.build();
    }

    @NotNull
    private static List<KnownHotspot> convertToHotspots(@NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry,
            @NotNull ProteinAnnotationExtractor proteinAnnotationExtractor) {
        List<KnownHotspot> hotspots = Lists.newArrayList();
        for (Map.Entry<ViccEntry, ViccExtractionResult> entryResult : resultsPerEntry.entrySet()) {
            ViccEntry entry = entryResult.getKey();
            for (Map.Entry<Feature, List<VariantHotspot>> featureResult : entryResult.getValue().hotspotsPerFeature().entrySet()) {
                Feature feature = featureResult.getKey();
                for (VariantHotspot hotspot : featureResult.getValue()) {
                    hotspots.add(ImmutableKnownHotspot.builder()
                            .from(hotspot)
                            .addSources(toKnowledgebase(entry.source()))
                            .gene(feature.geneSymbol())
                            .transcript(entry.transcriptId())
                            .proteinAnnotation(proteinAnnotationExtractor.apply(feature.name()))
                            .build());
                }
            }
        }

        return HotspotFunctions.consolidate(hotspots);
    }

    @NotNull
    private static List<KnownCopyNumber> convertToKnownAmpsDels(@NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
        List<KnownCopyNumber> copyNumbers = Lists.newArrayList();
        for (Map.Entry<ViccEntry, ViccExtractionResult> entry : resultsPerEntry.entrySet()) {
            copyNumbers.addAll(entry.getValue().ampsDelsPerFeature().values());
        }

        return CopyNumberFunctions.consolidate(copyNumbers);
    }

    @NotNull
    private static List<KnownFusionPair> convertToKnownFusions(@NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
        List<KnownFusionPair> fusions = Lists.newArrayList();
        for (Map.Entry<ViccEntry, ViccExtractionResult> entry : resultsPerEntry.entrySet()) {
            fusions.addAll(entry.getValue().fusionsPerFeature().values());
        }

        return FusionFunctions.consolidate(fusions);
    }

    private static void addActionability(@NotNull ImmutableExtractionOutput.Builder outputBuilder,
            @NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
        List<ActionableHotspot> actionableHotspots = Lists.newArrayList();
        List<ActionableRange> actionableRanges = Lists.newArrayList();
        List<ActionableGene> actionableGenes = Lists.newArrayList();
        List<ActionableFusion> actionableFusions = Lists.newArrayList();
        List<ActionableSignature> actionableSignatures = Lists.newArrayList();

        for (Map.Entry<ViccEntry, ViccExtractionResult> entry : resultsPerEntry.entrySet()) {
            ViccExtractionResult result = entry.getValue();
            ActionableEvent event = result.actionableEvent();
            if (event != null) {
                actionableHotspots.addAll(extractActionableHotspots(event, result.hotspotsPerFeature().values()));
                actionableRanges.addAll(extractActionableRanges(event,
                        result.geneRangesPerFeature().values(),
                        result.geneRangesPerFeature().keySet()));
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
    private static List<ActionableHotspot> extractActionableHotspots(@NotNull ActionableEvent actionableEvent,
            @NotNull Iterable<List<VariantHotspot>> hotspotLists) {
        List<ActionableHotspot> actionableHotspots = Lists.newArrayList();
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
    private static List<ActionableRange> extractActionableRanges(@NotNull ActionableEvent actionableEvent,
            @NotNull Iterable<List<GeneRangeAnnotation>> geneRangeAnnotationLists, @NotNull Set<Feature> features) {
        List<ActionableRange> actionableRanges = Lists.newArrayList();
        for (List<GeneRangeAnnotation> rangeList : geneRangeAnnotationLists) {
            for (GeneRangeAnnotation range : rangeList) {
                actionableRanges.add(ImmutableActionableRange.builder()
                        .from(actionableEvent)
                        .gene(range.gene())
                        .chromosome(range.chromosome())
                        .start(range.start())
                        .end(range.end())
                        .mutationType(range.mutationType())
                        .rangeInfo(range.rangeNumber())
                        .build());
            }
        }
        return actionableRanges;
    }

    @NotNull
    private static List<ActionableGene> extractActionableAmpsDels(@NotNull ActionableEvent actionableEvent,
            @NotNull Iterable<KnownCopyNumber> copyNumbers) {
        List<ActionableGene> actionableGenes = Lists.newArrayList();
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
    private static List<ActionableGene> extractActionableGeneLevelEvents(@NotNull ActionableEvent actionableEvent,
            @NotNull Iterable<GeneLevelAnnotation> geneLevelEvents) {
        List<ActionableGene> actionableGenes = Lists.newArrayList();
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
    private static List<ActionableFusion> extractActionableFusions(@NotNull ActionableEvent actionableEvent,
            @NotNull Iterable<KnownFusionPair> fusionAnnotations) {
        List<ActionableFusion> actionableFusions = Lists.newArrayList();
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
    private static List<ActionableSignature> extractActionableSignatures(@NotNull ActionableEvent actionableEvent,
            @NotNull Iterable<SignatureName> signatures) {
        List<ActionableSignature> actionableSignatures = Lists.newArrayList();
        for (SignatureName signature : signatures) {
            actionableSignatures.add(ImmutableActionableSignature.builder().from(actionableEvent).name(signature).build());
        }
        return actionableSignatures;
    }

    private static void printResults(@NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
        List<Feature> featuresWithoutGenomicEvents = Lists.newArrayList();
        int totalFeatureCount = 0;
        int featuresWithHotspotsCount = 0;
        int totalHotspotsCount = 0;
        int featuresWithCopyNumberCount = 0;
        int featuresWithFusionCount = 0;
        int featuresWithGeneLevelEventCount = 0;
        int featuresWithGeneRangeCount = 0;
        int totalRangeCount = 0;
        int featuresWithSignatureCount = 0;

        for (Map.Entry<ViccEntry, ViccExtractionResult> entry : resultsPerEntry.entrySet()) {
            ViccEntry viccEntry = entry.getKey();
            ViccExtractionResult viccExtractionResult = entry.getValue();
            for (Feature feature : viccEntry.features()) {
                List<VariantHotspot> hotspotsForFeature = viccExtractionResult.hotspotsPerFeature().get(feature);
                KnownCopyNumber ampDelForFeature = viccExtractionResult.ampsDelsPerFeature().get(feature);
                KnownFusionPair fusionForFeature = viccExtractionResult.fusionsPerFeature().get(feature);
                GeneLevelAnnotation geneLevelEventForFeature = viccExtractionResult.geneLevelEventsPerFeature().get(feature);
                List<GeneRangeAnnotation> geneRangesForFeature = viccExtractionResult.geneRangesPerFeature().get(feature);
                SignatureName signatureForFeature = viccExtractionResult.signaturesPerFeature().get(feature);

                if (hotspotsForFeature == null && ampDelForFeature == null && fusionForFeature == null && geneLevelEventForFeature == null
                        && geneRangesForFeature == null && signatureForFeature == null) {
                    if (feature.type() != MutationType.COMBINED && feature.type() != MutationType.COMPLEX) {
                        // For both combined and complex events we expect no genomic events to be derived.
                        featuresWithoutGenomicEvents.add(feature);
                    }
                } else {
                    if (hotspotsForFeature != null) {
                        featuresWithHotspotsCount++;
                        totalHotspotsCount += hotspotsForFeature.size();
                    }

                    if (ampDelForFeature != null) {
                        featuresWithCopyNumberCount++;
                    }

                    if (fusionForFeature != null) {
                        featuresWithFusionCount++;
                    }

                    if (geneLevelEventForFeature != null) {
                        featuresWithGeneLevelEventCount++;
                    }

                    if (geneRangesForFeature != null) {
                        featuresWithGeneRangeCount++;
                        totalRangeCount += geneRangesForFeature.size();
                    }

                    if (signatureForFeature != null) {
                        featuresWithSignatureCount++;
                    }
                }

                totalFeatureCount++;

            }
        }

        if (!featuresWithoutGenomicEvents.isEmpty()) {
            LOGGER.warn("No genomic events derived for {} features!", featuresWithoutGenomicEvents.size());
            for (Feature feature : featuresWithoutGenomicEvents) {
                LOGGER.debug(" No genomic events derived from '{}' in '{}'", feature.name(), feature.geneSymbol());
            }
        }

        LOGGER.info("Extraction performed on {} features in {} VICC entries", totalFeatureCount, resultsPerEntry.size());
        LOGGER.info(" Extracted {} hotspots for {} features", totalHotspotsCount, featuresWithHotspotsCount);
        LOGGER.info(" Extracted {} known amps and dels", featuresWithCopyNumberCount);
        LOGGER.info(" Extracted {} known fusions pairs", featuresWithFusionCount);
        LOGGER.info(" Extracted {} gene level events", featuresWithGeneLevelEventCount);
        LOGGER.info(" Extracted {} gene ranges for {} features", totalRangeCount, featuresWithGeneRangeCount);
        LOGGER.info(" Extracted {} signatures", featuresWithSignatureCount);
    }

    private static void writeInterpretationToTsv(@NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry,
            @NotNull String viccFeatureInterpretationTsv) throws IOException {
        List<String> lines = Lists.newArrayList();
        String header = new StringJoiner("\t").add("source").add("gene").add("event").add("type").add("interpretation").toString();
        lines.add(header);

        for (Map.Entry<ViccEntry, ViccExtractionResult> entry : resultsPerEntry.entrySet()) {
            ViccEntry viccEntry = entry.getKey();
            ViccExtractionResult viccExtractionResult = entry.getValue();
            for (Feature feature : viccEntry.features()) {
                StringJoiner interpretation = new StringJoiner(",");

                List<VariantHotspot> hotspotsForFeature = viccExtractionResult.hotspotsPerFeature().get(feature);
                if (hotspotsForFeature != null) {
                    interpretation.add(hotspotsForFeature.toString());
                }

                KnownCopyNumber ampDelForFeature = viccExtractionResult.ampsDelsPerFeature().get(feature);
                if (ampDelForFeature != null) {
                    interpretation.add(ampDelForFeature.toString());
                }

                KnownFusionPair fusionForFeature = viccExtractionResult.fusionsPerFeature().get(feature);
                if (fusionForFeature != null) {
                    interpretation.add(fusionForFeature.toString());
                }

                GeneLevelAnnotation geneLevelEventForFeature = viccExtractionResult.geneLevelEventsPerFeature().get(feature);
                if (geneLevelEventForFeature != null) {
                    interpretation.add(geneLevelEventForFeature.toString());
                }

                List<GeneRangeAnnotation> geneRangesForFeature = viccExtractionResult.geneRangesPerFeature().get(feature);
                if (geneRangesForFeature != null) {
                    interpretation.add(geneRangesForFeature.toString());
                }

                SignatureName signatureForFeature = viccExtractionResult.signaturesPerFeature().get(feature);
                if (signatureForFeature != null) {
                    interpretation.add(signatureForFeature.toString());
                }

                lines.add(new StringJoiner("\t").add(viccEntry.source().toString())
                        .add(feature.geneSymbol())
                        .add(feature.name())
                        .add(feature.type().name())
                        .add(interpretation.toString())
                        .toString());
            }

        }
        Files.write(new File(viccFeatureInterpretationTsv).toPath(), lines);
    }
}
