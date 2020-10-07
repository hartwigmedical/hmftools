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
import com.hartwig.hmftools.serve.hotspot.HotspotAnnotation;
import com.hartwig.hmftools.serve.hotspot.HotspotFunctions;
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
        List<ActionableFusion> actionableFusions = Lists.newArrayList();
        for (Map.Entry<ViccEntry, ViccExtractionResult> entry : resultsPerEntry.entrySet()) {
            Source source = fromViccSource(entry.getKey().source());
            ViccExtractionResult result = entry.getValue();
            ActionableEvidence evidence = result.actionableEvidence();
            if (evidence != null && source != null) {
                for (FusionAnnotation fusion : result.fusionsPerFeature().values()) {
                    if (fusion.fusionEvent() == FusionEvent.FUSION_PAIR) {
                        actionableFusions.add(ImmutableActionableFusion.builder()
                                .fusion(fusion.fusion())
                                .source(source)
                                .treatment(evidence.drugs())
                                .cancerType(evidence.cancerType())
                                .doid(evidence.doid())
                                .direction(evidence.direction())
                                .level(evidence.level())
                                .build());
                    }
                }
            }
        }
        String actionableFusionTsv = ActionableFusionFile.actionableFusionTsvPath(outputDir);
        LOGGER.info("Writing {} actionable fusions to {}", actionableFusions.size(), actionableFusionTsv);
        ActionableFusionFile.write(actionableFusionTsv, actionableFusions);
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
                String signatureForFeature = viccExtractionResult.signaturesPerFeature().get(feature);

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
                    events.add(signatureForFeature);
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
                String signatureForFeature = viccExtractionResult.signaturesPerFeature().get(feature);

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
