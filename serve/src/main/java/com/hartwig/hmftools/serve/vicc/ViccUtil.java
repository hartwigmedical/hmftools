package com.hartwig.hmftools.serve.vicc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotComparator;
import com.hartwig.hmftools.serve.Source;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusionFile;
import com.hartwig.hmftools.serve.actionability.fusion.ImmutableActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGeneFile;
import com.hartwig.hmftools.serve.actionability.gene.GeneEvent;
import com.hartwig.hmftools.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspotFile;
import com.hartwig.hmftools.serve.actionability.hotspot.ImmutableActionableHotspot;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.actionability.range.ActionableRangeFile;
import com.hartwig.hmftools.serve.actionability.range.ImmutableActionableRange;
import com.hartwig.hmftools.serve.actionability.signature.ActionableSignature;
import com.hartwig.hmftools.serve.actionability.signature.ActionableSignatureFile;
import com.hartwig.hmftools.serve.actionability.signature.ImmutableActionableSignature;
import com.hartwig.hmftools.serve.actionability.signature.SignatureName;
import com.hartwig.hmftools.serve.hotspot.HotspotAnnotation;
import com.hartwig.hmftools.serve.hotspot.HotspotFunctions;
import com.hartwig.hmftools.serve.vicc.copynumber.CopyNumberType;
import com.hartwig.hmftools.serve.vicc.copynumber.KnownAmplificationDeletion;
import com.hartwig.hmftools.serve.vicc.fusion.FusionAnnotation;
import com.hartwig.hmftools.serve.vicc.range.GeneRangeAnnotation;
import com.hartwig.hmftools.vicc.annotation.FusionEvent;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ViccUtil {

    private static final Logger LOGGER = LogManager.getLogger(ViccUtil.class);

    private static final String FIELD_DELIMITER = "\t";
    private static final String EVENT_DELIMITER = "|";

    private ViccUtil() {
    }

    public static void writeActionability(@NotNull String outputDir, @NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry)
            throws IOException {
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
                actionableGenes.addAll(extractActionablePromiscuousFusions(source, evidence, result.fusionsPerFeature().values()));
                actionableGenes.addAll(extractActionableAmpsDels(source, evidence, result.ampsDelsPerFeature().values()));
                actionableGenes.addAll(extractActionableGeneLevelEvents(source, evidence, result.geneLevelEventsPerFeature().values()));
                actionableFusions.addAll(extractActionableFusions(source, evidence, result.fusionsPerFeature().values()));
                actionableSignatures.addAll(extractActionableSignatures(source, evidence, result.signaturesPerFeature().values()));
            }
        }
        String actionableHotspotTsv = ActionableHotspotFile.actionableHotspotTsvPath(outputDir);
        LOGGER.info("Writing {} actionable hotspots to {}", actionableHotspots.size(), actionableHotspotTsv);
        ActionableHotspotFile.write(actionableHotspotTsv, actionableHotspots);

        String actionableRangeTsv = ActionableRangeFile.actionableRangeTsvPath(outputDir);
        LOGGER.info("Writing {} actionable ranges to {}", actionableRanges.size(), actionableRangeTsv);
        ActionableRangeFile.write(actionableRangeTsv, actionableRanges);

        String actionableGeneTsv = ActionableGeneFile.actionableGeneTsvPath(outputDir);
        LOGGER.info("Writing {} actionable genes to {}", actionableGenes.size(), actionableGeneTsv);
        ActionableGeneFile.write(actionableGeneTsv, actionableGenes);

        String actionableFusionTsv = ActionableFusionFile.actionableFusionTsvPath(outputDir);
        LOGGER.info("Writing {} actionable fusions to {}", actionableFusions.size(), actionableFusionTsv);
        ActionableFusionFile.write(actionableFusionTsv, actionableFusions);

        String actionableSignatureTsv = ActionableSignatureFile.actionableSignatureTsvPath(outputDir);
        LOGGER.info("Writing {} actionable signatures to {}", actionableSignatures.size(), actionableSignatureTsv);
        ActionableSignatureFile.write(actionableSignatureTsv, actionableSignatures);
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
                        .build());
            }
        }
        return actionableRanges;
    }

    @NotNull
    private static List<ActionableGene> extractActionablePromiscuousFusions(@NotNull Source source, @NotNull ActionableEvidence evidence,
            @NotNull Iterable<FusionAnnotation> fusionAnnotations) {
        List<ActionableGene> actionableGenes = Lists.newArrayList();
        for (FusionAnnotation fusion : fusionAnnotations) {
            if (fusion.fusionEvent() == FusionEvent.FUSION_PROMISCUOUS) {
                actionableGenes.add(ImmutableActionableGene.builder()
                        .gene(fusion.fusion())
                        .event(GeneEvent.FUSION)
                        .source(source)
                        .treatment(evidence.drugs())
                        .cancerType(evidence.cancerType())
                        .doid(evidence.doid())
                        .direction(evidence.direction())
                        .level(evidence.level())
                        .build());
            }
        }
        return actionableGenes;
    }

    @NotNull
    private static List<ActionableGene> extractActionableAmpsDels(@NotNull Source source, @NotNull ActionableEvidence evidence,
            @NotNull Iterable<KnownAmplificationDeletion> ampsDels) {
        List<ActionableGene> actionableGenes = Lists.newArrayList();
        for (KnownAmplificationDeletion ampDel : ampsDels) {
            actionableGenes.add(ImmutableActionableGene.builder()
                    // TODO Improve event mapping
                    .gene(ampDel.gene())
                    .event(ampDel.type() == CopyNumberType.AMPLIFICATION ? GeneEvent.AMPLIFICATION : GeneEvent.DELETION)
                    .source(source)
                    .treatment(evidence.drugs())
                    .cancerType(evidence.cancerType())
                    .doid(evidence.doid())
                    .direction(evidence.direction())
                    .level(evidence.level())
                    .build());
        }
        return actionableGenes;
    }

    @NotNull
    private static List<ActionableGene> extractActionableGeneLevelEvents(@NotNull Source source, @NotNull ActionableEvidence evidence,
            @NotNull Iterable<String> geneLevelEvents) {
        List<ActionableGene> actionableGenes = Lists.newArrayList();
        for (String geneLevelEvent : geneLevelEvents) {
            actionableGenes.add(ImmutableActionableGene.builder()
                    // TODO Implement event
                    .gene(geneLevelEvent)
                    .event(GeneEvent.ACTIVATION)
                    .source(source)
                    .treatment(evidence.drugs())
                    .cancerType(evidence.cancerType())
                    .doid(evidence.doid())
                    .direction(evidence.direction())
                    .level(evidence.level())
                    .build());
        }
        return actionableGenes;
    }

    @NotNull
    private static List<ActionableFusion> extractActionableFusions(@NotNull Source source, @NotNull ActionableEvidence evidence,
            @NotNull Iterable<FusionAnnotation> fusionAnnotations) {
        List<ActionableFusion> actionableFusions = Lists.newArrayList();
        for (FusionAnnotation fusion : fusionAnnotations) {
            if (fusion.fusionEvent() == FusionEvent.FUSION_PAIR) {
                actionableFusions.add(ImmutableActionableFusion.builder()
                        // TODO Separate gene up from gene down.
                        .geneUp(fusion.fusion())
                        .exonUp(fusion.exonUp())
                        .geneDown(fusion.fusion())
                        .exonDown(fusion.exonDown())
                        .source(source)
                        .treatment(evidence.drugs())
                        .cancerType(evidence.cancerType())
                        .doid(evidence.doid())
                        .direction(evidence.direction())
                        .level(evidence.level())
                        .build());
            }
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

    public static void writeFeatureTypes(@NotNull String featureTypeTsv, @NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry)
            throws IOException {
        List<String> lines = Lists.newArrayList();
        String header = new StringJoiner(FIELD_DELIMITER).add("name").add("gene").add("type").add("biomarkerType").add("events").toString();
        lines.add(header);

        for (Map.Entry<ViccEntry, ViccExtractionResult> entry : resultsPerEntry.entrySet()) {
            ViccEntry viccEntry = entry.getKey();
            ViccExtractionResult viccExtractionResult = entry.getValue();
            for (Feature feature : viccEntry.features()) {
                List<VariantHotspot> hotspotsForFeature = viccExtractionResult.hotspotsPerFeature().get(feature);
                KnownAmplificationDeletion ampDelForFeature = viccExtractionResult.ampsDelsPerFeature().get(feature);
                FusionAnnotation fusionForFeature = viccExtractionResult.fusionsPerFeature().get(feature);
                String geneLevelEventForFeature = viccExtractionResult.geneLevelEventsPerFeature().get(feature);
                List<GeneRangeAnnotation> geneRangeForFeature = viccExtractionResult.geneRangesPerFeature().get(feature);
                SignatureName signatureForFeature = viccExtractionResult.signaturesPerFeature().get(feature);

                StringJoiner events = new StringJoiner(EVENT_DELIMITER);
                if (hotspotsForFeature != null) {
                    for (VariantHotspot hotspot : hotspotsForFeature) {
                        events.add(hotspot.toString());
                    }
                }

                if (ampDelForFeature != null) {
                    events.add(ampDelForFeature.toString());
                }

                if (fusionForFeature != null) {
                    events.add(fusionForFeature.toString());
                }

                if (geneLevelEventForFeature != null) {
                    events.add(geneLevelEventForFeature);
                }

                if (geneRangeForFeature != null) {
                    events.add(geneRangeForFeature.toString());
                }

                if (signatureForFeature != null) {
                    events.add(signatureForFeature.toString());
                }

                StringJoiner featureString = new StringJoiner(FIELD_DELIMITER);
                featureString.add(feature.name())
                        .add(feature.geneSymbol())
                        .add(feature.type().toString())
                        .add(feature.biomarkerType())
                        .add(events.toString());
                lines.add(featureString.toString());
            }
        }

        LOGGER.info("Writing {} feature types to {}", lines.size() - 1, featureTypeTsv);
        Files.write(new File(featureTypeTsv).toPath(), lines);
    }

    @NotNull
    public static Map<VariantHotspot, HotspotAnnotation> convertToHotspotMap(
            @NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
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

    public static void printResults(@NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
        List<Feature> featuresWithoutGenomicEvents = Lists.newArrayList();
        int totalFeatureCount = 0;
        int featuresWithHotspotsCount = 0;
        int totalHotspotsCount = 0;
        int featuresWithCopyNumberCount = 0;
        int featuresWithFusionCount = 0;
        int featuresWithGeneLevelEventCount = 0;
        int featuresWithGeneRangeCount = 0;
        int featuresWithSignatureCount = 0;

        for (Map.Entry<ViccEntry, ViccExtractionResult> entry : resultsPerEntry.entrySet()) {
            ViccEntry viccEntry = entry.getKey();
            ViccExtractionResult viccExtractionResult = entry.getValue();
            for (Feature feature : viccEntry.features()) {
                List<VariantHotspot> hotspotsForFeature = viccExtractionResult.hotspotsPerFeature().get(feature);
                KnownAmplificationDeletion ampDelForFeature = viccExtractionResult.ampsDelsPerFeature().get(feature);
                FusionAnnotation fusionForFeature = viccExtractionResult.fusionsPerFeature().get(feature);
                String geneLevelEventForFeature = viccExtractionResult.geneLevelEventsPerFeature().get(feature);
                List<GeneRangeAnnotation> geneRangeForFeature = viccExtractionResult.geneRangesPerFeature().get(feature);
                SignatureName signatureForFeature = viccExtractionResult.signaturesPerFeature().get(feature);

                if (hotspotsForFeature == null && ampDelForFeature == null && fusionForFeature == null && geneLevelEventForFeature == null
                        && geneRangeForFeature == null && signatureForFeature == null) {
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

                    if (geneRangeForFeature != null) {
                        featuresWithGeneRangeCount++;
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
        LOGGER.info(" Extracted {} gene ranges", featuresWithGeneRangeCount);
        LOGGER.info(" Extracted {} signatures", featuresWithSignatureCount);
    }
}
