package com.hartwig.hmftools.serve.vicc;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotComparator;
import com.hartwig.hmftools.serve.Source;
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
import com.hartwig.hmftools.serve.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.hotspot.HotspotAnnotation;
import com.hartwig.hmftools.serve.hotspot.HotspotFunctions;
import com.hartwig.hmftools.serve.vicc.copynumber.CopyNumberExtractor;
import com.hartwig.hmftools.serve.vicc.fusion.FusionExtractor;
import com.hartwig.hmftools.serve.vicc.genelevel.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.vicc.genelevel.GeneLevelEventExtractor;
import com.hartwig.hmftools.serve.vicc.hotspot.HotspotExtractor;
import com.hartwig.hmftools.serve.vicc.range.GeneRangeAnnotation;
import com.hartwig.hmftools.serve.vicc.range.GeneRangeExtractor;
import com.hartwig.hmftools.serve.vicc.signatures.SignaturesExtractor;
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
    private final HotspotExtractor hotspotExtractor;
    @NotNull
    private final CopyNumberExtractor copyNumberExtractor;
    @NotNull
    private final FusionExtractor fusionExtractor;
    @NotNull
    private final GeneLevelEventExtractor geneLevelEventExtractor;
    @NotNull
    private final GeneRangeExtractor geneRangeExtractor;
    @NotNull
    private final SignaturesExtractor signaturesExtractor;

    public ViccExtractor(@NotNull final HotspotExtractor hotspotExtractor, @NotNull final CopyNumberExtractor copyNumberExtractor,
            @NotNull final FusionExtractor fusionExtractor, @NotNull final GeneLevelEventExtractor geneLevelEventExtractor,
            @NotNull final GeneRangeExtractor geneRangeExtractor, @NotNull final SignaturesExtractor signaturesExtractor) {
        this.hotspotExtractor = hotspotExtractor;
        this.copyNumberExtractor = copyNumberExtractor;
        this.fusionExtractor = fusionExtractor;
        this.geneLevelEventExtractor = geneLevelEventExtractor;
        this.geneRangeExtractor = geneRangeExtractor;
        this.signaturesExtractor = signaturesExtractor;
    }

    @NotNull
    public ViccExtractionOutput extractFromViccEntries(@NotNull List<ViccEntry> viccEntries) {
        Map<ViccEntry, ViccExtractionResult> resultsPerEntry = Maps.newHashMap();
        for (ViccEntry entry : viccEntries) {
            Map<Feature, List<VariantHotspot>> hotspotsPerFeature = hotspotExtractor.extractHotspots(entry);
            Map<Feature, KnownCopyNumber> ampsDelsPerFeature = copyNumberExtractor.extractKnownAmplificationsDeletions(entry);
            Map<Feature, KnownFusionPair> fusionsPerFeature = fusionExtractor.extractKnownFusions(entry);
            Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = geneLevelEventExtractor.extractKnownGeneLevelEvents(entry);
            Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature = geneRangeExtractor.extractGeneRanges(entry);
            Map<Feature, SignatureName> signaturesPerFeature = signaturesExtractor.extractSignatures(entry);

            ActionableEvidence actionableEvidence = ActionableEvidenceFactory.toActionableEvidence(entry);

            resultsPerEntry.put(entry,
                    ImmutableViccExtractionResult.builder()
                            .hotspotsPerFeature(hotspotsPerFeature)
                            .ampsDelsPerFeature(ampsDelsPerFeature)
                            .fusionsPerFeature(fusionsPerFeature)
                            .geneLevelEventsPerFeature(geneLevelEventsPerFeature)
                            .geneRangesPerFeature(geneRangesPerFeature)
                            .signaturesPerFeature(signaturesPerFeature)
                            .actionableEvidence(actionableEvidence)
                            .build());
        }
        printResults(resultsPerEntry);

        ImmutableViccExtractionOutput.Builder outputBuilder = ImmutableViccExtractionOutput.builder()
                .hotspots(convertToHotspots(resultsPerEntry))
                .knownAmpsDels(convertToKnownAmpsDels(resultsPerEntry))
                .knownFusions(convertToKnownFusions(resultsPerEntry));

        addActionability(outputBuilder, resultsPerEntry);

        return outputBuilder.build();
    }

    @NotNull
    private static Map<VariantHotspot, HotspotAnnotation> convertToHotspots(@NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
        Map<VariantHotspot, HotspotAnnotation> convertedMap = Maps.newTreeMap(new VariantHotspotComparator());
        for (Map.Entry<ViccEntry, ViccExtractionResult> entryResult : resultsPerEntry.entrySet()) {
            ViccEntry entry = entryResult.getKey();
            for (Map.Entry<Feature, List<VariantHotspot>> featureResult : entryResult.getValue().hotspotsPerFeature().entrySet()) {
                Feature feature = featureResult.getKey();
                for (VariantHotspot hotspot : featureResult.getValue()) {
                    HotspotAnnotation currentAnnotation = convertedMap.get(hotspot);
                    HotspotAnnotation newAnnotation = new HotspotAnnotation(Sets.newHashSet(entry.source().display()),
                            feature.geneSymbol(),
                            entry.transcriptId(),
                            feature.proteinAnnotation());
                    if (currentAnnotation != null) {
                        convertedMap.put(hotspot, HotspotFunctions.mergeHotspotAnnotations(currentAnnotation, newAnnotation));
                    } else {
                        convertedMap.put(hotspot, newAnnotation);
                    }

                }
            }
        }

        return convertedMap;
    }

    @NotNull
    private static List<KnownCopyNumber> convertToKnownAmpsDels(@NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
        List<KnownCopyNumber> copyNumbers = Lists.newArrayList();
        for (Map.Entry<ViccEntry, ViccExtractionResult> entry : resultsPerEntry.entrySet()) {
            copyNumbers.addAll(entry.getValue().ampsDelsPerFeature().values());
        }
        return copyNumbers;
    }

    @NotNull
    private static List<KnownFusionPair> convertToKnownFusions(@NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
        List<KnownFusionPair> fusions = Lists.newArrayList();
        for (Map.Entry<ViccEntry, ViccExtractionResult> entry : resultsPerEntry.entrySet()) {
            fusions.addAll(entry.getValue().fusionsPerFeature().values());
        }
        return fusions;
    }

    private static void addActionability(@NotNull ImmutableViccExtractionOutput.Builder outputBuilder,
            @NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
        List<ActionableHotspot> actionableHotspots = Lists.newArrayList();
        List<ActionableRange> actionableRanges = Lists.newArrayList();
        List<ActionableGene> actionableGenes = Lists.newArrayList();
        List<ActionableFusion> actionableFusions = Lists.newArrayList();
        List<ActionableSignature> actionableSignatures = Lists.newArrayList();

        for (Map.Entry<ViccEntry, ViccExtractionResult> entry : resultsPerEntry.entrySet()) {
            Source source = fromViccSource(entry.getKey().source());
            ViccExtractionResult result = entry.getValue();
            ActionableEvidence evidence = result.actionableEvidence();
            if (evidence != null && source != null) {
                actionableHotspots.addAll(extractActionableHotspots(source, evidence, result.hotspotsPerFeature().values()));
                actionableRanges.addAll(extractActionableRanges(source, evidence, result.geneRangesPerFeature().values()));
                actionableGenes.addAll(extractActionableAmpsDels(source, evidence, result.ampsDelsPerFeature().values()));
                actionableGenes.addAll(extractActionableGeneLevelEvents(source, evidence, result.geneLevelEventsPerFeature().values()));
                actionableFusions.addAll(extractActionableFusions(source, evidence, result.fusionsPerFeature().values()));
                actionableSignatures.addAll(extractActionableSignatures(source, evidence, result.signaturesPerFeature().values()));
            }
        }

        outputBuilder.actionableHotspots(actionableHotspots);
        outputBuilder.actionableRanges(actionableRanges);
        outputBuilder.actionableGenes(actionableGenes);
        outputBuilder.actionableFusions(actionableFusions);
        outputBuilder.actionableSignatures(actionableSignatures);
    }

    @NotNull
    private static List<ActionableHotspot> extractActionableHotspots(@NotNull Source source, @NotNull ActionableEvidence evidence,
            @NotNull Iterable<List<VariantHotspot>> hotspotLists) {
        List<ActionableHotspot> actionableHotspots = Lists.newArrayList();
        for (List<VariantHotspot> hotspotList : hotspotLists) {
            for (VariantHotspot hotspot : hotspotList) {
                actionableHotspots.add(ImmutableActionableHotspot.builder()
                        .chromosome(hotspot.chromosome())
                        .position(hotspot.position())
                        .ref(hotspot.ref())
                        .alt(hotspot.alt())
                        .source(source)
                        .treatment(evidence.drugs())
                        .cancerType(evidence.cancerType())
                        .doid(evidence.doid())
                        .direction(evidence.direction())
                        .level(evidence.level())
                        .url(evidence.url())
                        .build());
            }
        }
        return actionableHotspots;
    }

    @NotNull
    private static List<ActionableRange> extractActionableRanges(@NotNull Source source, @NotNull ActionableEvidence evidence,
            @NotNull Iterable<List<GeneRangeAnnotation>> geneRangeAnnotationLists) {
        List<ActionableRange> actionableRanges = Lists.newArrayList();
        for (List<GeneRangeAnnotation> rangeList : geneRangeAnnotationLists) {
            for (GeneRangeAnnotation range : rangeList) {
                actionableRanges.add(ImmutableActionableRange.builder()
                        .gene(range.gene())
                        .chromosome(range.chromosome())
                        .start(range.start())
                        .end(range.end())
                        .mutationType(range.mutationType())
                        .source(source)
                        .treatment(evidence.drugs())
                        .cancerType(evidence.cancerType())
                        .doid(evidence.doid())
                        .direction(evidence.direction())
                        .level(evidence.level())
                        .url(evidence.url())
                        .build());
            }
        }
        return actionableRanges;
    }

    @NotNull
    private static List<ActionableGene> extractActionableAmpsDels(@NotNull Source source, @NotNull ActionableEvidence evidence,
            @NotNull Iterable<KnownCopyNumber> ampsDels) {
        List<ActionableGene> actionableGenes = Lists.newArrayList();
        for (KnownCopyNumber ampDel : ampsDels) {
            GeneLevelEvent event;
            switch (ampDel.type()) {
                case AMPLIFICATION: {
                    event = GeneLevelEvent.AMPLIFICATION;
                    break;
                }
                case DELETION: {
                    event = GeneLevelEvent.DELETION;
                    break;
                }
                default:
                    throw new IllegalStateException("Invalid copy number type: " + ampDel.type());
            }

            actionableGenes.add(ImmutableActionableGene.builder()
                    .gene(ampDel.gene())
                    .event(event)
                    .source(source)
                    .treatment(evidence.drugs())
                    .cancerType(evidence.cancerType())
                    .doid(evidence.doid())
                    .direction(evidence.direction())
                    .level(evidence.level())
                    .url(evidence.url())
                    .build());
        }
        return actionableGenes;
    }

    @NotNull
    private static List<ActionableGene> extractActionableGeneLevelEvents(@NotNull Source source, @NotNull ActionableEvidence evidence,
            @NotNull Iterable<GeneLevelAnnotation> geneLevelEvents) {
        List<ActionableGene> actionableGenes = Lists.newArrayList();
        for (GeneLevelAnnotation geneLevelEvent : geneLevelEvents) {
            actionableGenes.add(ImmutableActionableGene.builder()
                    .gene(geneLevelEvent.gene())
                    .event(geneLevelEvent.event())
                    .source(source)
                    .treatment(evidence.drugs())
                    .cancerType(evidence.cancerType())
                    .doid(evidence.doid())
                    .direction(evidence.direction())
                    .level(evidence.level())
                    .url(evidence.url())
                    .build());
        }
        return actionableGenes;
    }

    @NotNull
    private static List<ActionableFusion> extractActionableFusions(@NotNull Source source, @NotNull ActionableEvidence evidence,
            @NotNull Iterable<KnownFusionPair> fusionAnnotations) {
        List<ActionableFusion> actionableFusions = Lists.newArrayList();
        for (KnownFusionPair fusion : fusionAnnotations) {
            actionableFusions.add(ImmutableActionableFusion.builder()
                    .geneUp(fusion.geneUp())
                    .exonUp(fusion.exonUp())
                    .geneDown(fusion.geneDown())
                    .exonDown(fusion.exonDown())
                    .source(source)
                    .treatment(evidence.drugs())
                    .cancerType(evidence.cancerType())
                    .doid(evidence.doid())
                    .direction(evidence.direction())
                    .level(evidence.level())
                    .url(evidence.url())
                    .build());
        }
        return actionableFusions;
    }

    @NotNull
    private static List<ActionableSignature> extractActionableSignatures(@NotNull Source source, @NotNull ActionableEvidence evidence,
            @NotNull Iterable<SignatureName> signatures) {
        List<ActionableSignature> actionableSignatures = Lists.newArrayList();
        for (SignatureName signature : signatures) {
            actionableSignatures.add(ImmutableActionableSignature.builder()
                    .name(signature)
                    .source(source)
                    .treatment(evidence.drugs())
                    .cancerType(evidence.cancerType())
                    .doid(evidence.doid())
                    .direction(evidence.direction())
                    .level(evidence.level())
                    .build());
        }
        return actionableSignatures;
    }

    @Nullable
    private static Source fromViccSource(@NotNull ViccSource source) {
        switch (source) {
            case CIVIC:
                return Source.CIVIC;
            case CGI:
                return Source.CGI;
            case JAX:
                return Source.JAX;
            case ONCOKB:
                return Source.ONCOKB;
            default:
                return null;
        }
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
                    featuresWithoutGenomicEvents.add(feature);
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

        LOGGER.info("No genomic events derived for {} features.", featuresWithoutGenomicEvents.size());
        for (Feature feature : featuresWithoutGenomicEvents) {
            LOGGER.debug(" No genomic events derived from '{}' in '{}'", feature.name(), feature.geneSymbol());
        }

        LOGGER.info("Extraction performed on {} features from {} entries", totalFeatureCount, resultsPerEntry.size());
        LOGGER.info(" Extracted {} hotspots for {} features", totalHotspotsCount, featuresWithHotspotsCount);
        LOGGER.info(" Extracted {} known amps and dels", featuresWithCopyNumberCount);
        LOGGER.info(" Extracted {} fusions", featuresWithFusionCount);
        LOGGER.info(" Extracted {} gene level events", featuresWithGeneLevelEventCount);
        LOGGER.info(" Extracted {} gene ranges for {} features", totalRangeCount, featuresWithGeneRangeCount);
        LOGGER.info(" Extracted {} signatures", featuresWithSignatureCount);
    }
}
