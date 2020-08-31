package com.hartwig.hmftools.serve.vicc;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotComparator;
import com.hartwig.hmftools.serve.hotspot.HotspotAnnotation;
import com.hartwig.hmftools.serve.hotspot.HotspotFunctions;
import com.hartwig.hmftools.serve.vicc.copynumber.KnownAmplificationDeletion;
import com.hartwig.hmftools.serve.vicc.fusion.FusionAnnotation;
import com.hartwig.hmftools.serve.vicc.range.GeneRangeAnnotation;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ViccUtil {

    private static final Logger LOGGER = LogManager.getLogger(ViccUtil.class);

    private ViccUtil() {
    }

    public static void writeEventMappingToTsv(@NotNull String eventMappingTsv,
            @NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(eventMappingTsv));
        writer.write("Mapped event" + "\t" + "Gene" + "\t" + "Name" + "\t" + "Event" + "\t" + "Feature" + "\n");

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
                        writer.write("Hotspot" + "\t" + feature.geneSymbol() + "\t" + feature.name() + "\t" + feature.biomarkerType() + "\t"
                                + feature + "\n");
                    }

                    if (ampDelForFeature != null) {
                        writer.write("Amp/Del" + "\t" + feature.geneSymbol() + "\t" + feature.name() + "\t" + feature.biomarkerType() + "\t"
                                + feature + "\n");

                        featuresWithCopyNumberCount++;
                    }

                    if (fusionForFeature != null) {
                        writer.write(fusionForFeature.fusionEvent() + "\t" + fusionForFeature.fusion() + "\t" + feature.name() + "\t"
                                + feature.biomarkerType() + "\t" + feature + "\n");
                        featuresWithFusionCount++;
                    }

                    if (geneLevelEventForFeature != null) {
                        writer.write(
                                "Gene_level" + "\t" + feature.geneSymbol() + "\t" + feature.name() + "\t" + feature.biomarkerType() + "\t"
                                        + feature + "\n");
                        featuresWithGeneLevelEventCount++;
                    }

                    if (geneRangeForFeature != null) {
                        writer.write(
                                "Gene_range" + "\t" + feature.geneSymbol() + "\t" + feature.name() + "\t" + feature.biomarkerType() + "\t"
                                        + feature + "\n");
                        featuresWithGeneRangeCount++;
                    }

                    if (signatureForFeature != null) {
                        writer.write(
                                "Signature" + "\t" + feature.geneSymbol() + "\t" + feature.name() + "\t" + feature.biomarkerType() + "\t"
                                        + feature + "\n");

                        featuresWithSignatureCount++;
                    }
                }

                totalFeatureCount++;
            }
        }

        LOGGER.info("No genomic events derived for {} features.", featuresWithoutGenomicEvents.size());
        for (Feature feature : featuresWithoutGenomicEvents) {
            if (!FeatureIgnoreUtil.canIgnore(feature)) {
                LOGGER.debug(" No genomic events derived from '{}' in '{}'", feature.name(), feature.geneSymbol());
                writer.write("Unmapped_event" + "\t" + feature.geneSymbol() + "\t" + feature.name() + "\t" + feature.biomarkerType() + "\t"
                        + feature + "\n");
            } else {
                writer.write("Ignore_event" + "\t" + feature.geneSymbol() + "\t" + feature.name() + "\t" + feature.biomarkerType() + "\t"
                        + feature + "\n");
            }
        }

        LOGGER.info("Extraction performed on {} features from {} entries", totalFeatureCount, resultsPerEntry.size());
        LOGGER.info(" Extracted {} hotspots for {} features", totalHotspotsCount, featuresWithHotspotsCount);
        LOGGER.info(" Extracted {} known amps and dels", featuresWithCopyNumberCount);
        LOGGER.info(" Extracted {} fusions", featuresWithFusionCount);
        LOGGER.info(" Extracted {} gene level events", featuresWithGeneLevelEventCount);
        LOGGER.info(" Extracted {} gene ranges", featuresWithGeneRangeCount);
        LOGGER.info(" Extracted {} signatures", featuresWithSignatureCount);

        writer.close();
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
}
