package com.hartwig.hmftools.serve.vicc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusionFile;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGeneFile;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspotFile;
import com.hartwig.hmftools.serve.actionability.range.ActionableRangeFile;
import com.hartwig.hmftools.serve.actionability.signature.ActionableSignatureFile;
import com.hartwig.hmftools.serve.actionability.signature.SignatureName;
import com.hartwig.hmftools.serve.vicc.copynumber.CopyNumberAnnotation;
import com.hartwig.hmftools.serve.vicc.fusion.FusionAnnotation;
import com.hartwig.hmftools.serve.vicc.genelevel.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.vicc.range.GeneRangeAnnotation;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ViccUtil {

    private static final Logger LOGGER = LogManager.getLogger(ViccUtil.class);

    private static final String FIELD_DELIMITER = "\t";
    private static final String EVENT_DELIMITER = "|";

    private ViccUtil() {
    }

    public static void writeActionability(@NotNull String outputDir, @NotNull ViccExtractionOutput viccExtractionOutput)
            throws IOException {
        String actionableHotspotTsv = ActionableHotspotFile.actionableHotspotTsvPath(outputDir);
        LOGGER.info("Writing {} actionable hotspots to {}", viccExtractionOutput.actionableHotspots().size(), actionableHotspotTsv);
        ActionableHotspotFile.write(actionableHotspotTsv, viccExtractionOutput.actionableHotspots());

        String actionableRangeTsv = ActionableRangeFile.actionableRangeTsvPath(outputDir);
        LOGGER.info("Writing {} actionable ranges to {}", viccExtractionOutput.actionableRanges().size(), actionableRangeTsv);
        ActionableRangeFile.write(actionableRangeTsv, viccExtractionOutput.actionableRanges());

        String actionableGeneTsv = ActionableGeneFile.actionableGeneTsvPath(outputDir);
        LOGGER.info("Writing {} actionable genes to {}", viccExtractionOutput.actionableGenes().size(), actionableGeneTsv);
        ActionableGeneFile.write(actionableGeneTsv, viccExtractionOutput.actionableGenes());

        String actionableFusionTsv = ActionableFusionFile.actionableFusionTsvPath(outputDir);
        LOGGER.info("Writing {} actionable fusions to {}", viccExtractionOutput.actionableFusions().size(), actionableFusionTsv);
        ActionableFusionFile.write(actionableFusionTsv, viccExtractionOutput.actionableFusions());

        String actionableSignatureTsv = ActionableSignatureFile.actionableSignatureTsvPath(outputDir);
        LOGGER.info("Writing {} actionable signatures to {}", viccExtractionOutput.actionableSignatures().size(), actionableSignatureTsv);
        ActionableSignatureFile.write(actionableSignatureTsv, viccExtractionOutput.actionableSignatures());
    }

    public static void writeFeatureTypes(@NotNull String featureTypeTsv, @NotNull ViccExtractionOutput viccExtractionOutput)
            throws IOException {
        List<String> lines = Lists.newArrayList();
        String header = new StringJoiner(FIELD_DELIMITER).add("name").add("gene").add("type").add("biomarkerType").add("events").toString();
        lines.add(header);

        for (Map.Entry<ViccEntry, ViccExtractionResult> entry : viccExtractionOutput.resultsPerEntry().entrySet()) {
            ViccEntry viccEntry = entry.getKey();
            ViccExtractionResult viccExtractionResult = entry.getValue();
            for (Feature feature : viccEntry.features()) {
                List<VariantHotspot> hotspotsForFeature = viccExtractionResult.hotspotsPerFeature().get(feature);
                CopyNumberAnnotation ampDelForFeature = viccExtractionResult.ampsDelsPerFeature().get(feature);
                FusionAnnotation fusionForFeature = viccExtractionResult.fusionsPerFeature().get(feature);
                GeneLevelAnnotation geneLevelEventForFeature = viccExtractionResult.geneLevelEventsPerFeature().get(feature);
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
                    events.add(geneLevelEventForFeature.toString());
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

    public static void printResults(@NotNull ViccExtractionOutput viccExtractionOutput) {
        List<Feature> featuresWithoutGenomicEvents = Lists.newArrayList();
        int totalFeatureCount = 0;
        int featuresWithHotspotsCount = 0;
        int totalHotspotsCount = 0;
        int featuresWithCopyNumberCount = 0;
        int featuresWithFusionCount = 0;
        int featuresWithGeneLevelEventCount = 0;
        int featuresWithGeneRangeCount = 0;
        int featuresWithSignatureCount = 0;

        for (Map.Entry<ViccEntry, ViccExtractionResult> entry : viccExtractionOutput.resultsPerEntry().entrySet()) {
            ViccEntry viccEntry = entry.getKey();
            ViccExtractionResult viccExtractionResult = entry.getValue();
            for (Feature feature : viccEntry.features()) {
                List<VariantHotspot> hotspotsForFeature = viccExtractionResult.hotspotsPerFeature().get(feature);
                CopyNumberAnnotation ampDelForFeature = viccExtractionResult.ampsDelsPerFeature().get(feature);
                FusionAnnotation fusionForFeature = viccExtractionResult.fusionsPerFeature().get(feature);
                GeneLevelAnnotation geneLevelEventForFeature = viccExtractionResult.geneLevelEventsPerFeature().get(feature);
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

        LOGGER.info("Extraction performed on {} features from {} entries",
                totalFeatureCount,
                viccExtractionOutput.resultsPerEntry().size());
        LOGGER.info(" Extracted {} hotspots for {} features", totalHotspotsCount, featuresWithHotspotsCount);
        LOGGER.info(" Extracted {} known amps and dels", featuresWithCopyNumberCount);
        LOGGER.info(" Extracted {} fusions", featuresWithFusionCount);
        LOGGER.info(" Extracted {} gene level events", featuresWithGeneLevelEventCount);
        LOGGER.info(" Extracted {} gene ranges", featuresWithGeneRangeCount);
        LOGGER.info(" Extracted {} signatures", featuresWithSignatureCount);
    }
}
