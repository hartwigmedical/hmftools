package com.hartwig.hmftools.serve.sources.vicc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.extraction.codon.CodonAnnotation;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.exon.ExonAnnotation;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.extraction.signature.SignatureName;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ViccUtil {

    private static final Logger LOGGER = LogManager.getLogger(ViccUtil.class);

    private static final String FIELD_DELIMITER = "\t";
    private static final String SUB_FIELD_DELIMITER = ",";

    private ViccUtil() {
    }

    public static void printExtractionResults(@NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
        List<Feature> featuresWithoutGenomicEvents = Lists.newArrayList();

        int totalFeatureCount = 0;
        int featuresWithHotspotsCount = 0;
        int totalHotspotsCount = 0;
        int featuresWithCodonCount = 0;
        int totalCodonCount = 0;
        int featuresWithExonCount = 0;
        int totalExonCount = 0;

        int featuresWithGeneLevelEventCount = 0;
        int featuresWithCopyNumberCount = 0;
        int featuresWithFusionCount = 0;
        int featuresWithSignatureCount = 0;

        for (Map.Entry<ViccEntry, ViccExtractionResult> resultPerEntry : resultsPerEntry.entrySet()) {
            ViccEntry entry = resultPerEntry.getKey();
            ViccExtractionResult result = resultPerEntry.getValue();
            for (Feature feature : entry.features()) {
                List<VariantHotspot> hotspotsForFeature = result.hotspotsPerFeature().get(feature);
                List<CodonAnnotation> codonsForFeature = result.codonsPerFeature().get(feature);
                List<ExonAnnotation> exonsForFeature = result.exonsPerFeature().get(feature);
                GeneLevelAnnotation geneLevelEventForFeature = result.geneLevelEventsPerFeature().get(feature);
                KnownCopyNumber ampDelForFeature = result.ampsDelsPerFeature().get(feature);
                KnownFusionPair fusionForFeature = result.fusionsPerFeature().get(feature);
                SignatureName signatureForFeature = result.signaturesPerFeature().get(feature);

                if (hotspotsForFeature == null && codonsForFeature == null && exonsForFeature == null && geneLevelEventForFeature == null
                        && ampDelForFeature == null && fusionForFeature == null && signatureForFeature == null) {
                    if (feature.type() != EventType.COMBINED && feature.type() != EventType.COMPLEX) {
                        // For both combined and complex events we expect no genomic events to be derived.
                        featuresWithoutGenomicEvents.add(feature);
                    }
                } else {
                    if (hotspotsForFeature != null) {
                        featuresWithHotspotsCount++;
                        totalHotspotsCount += hotspotsForFeature.size();
                    }

                    if (codonsForFeature != null) {
                        featuresWithCodonCount++;
                        totalCodonCount += codonsForFeature.size();
                    }

                    if (exonsForFeature != null) {
                        featuresWithExonCount++;
                        totalExonCount += exonsForFeature.size();
                    }

                    if (geneLevelEventForFeature != null) {
                        featuresWithGeneLevelEventCount++;
                    }

                    if (ampDelForFeature != null) {
                        featuresWithCopyNumberCount++;
                    }

                    if (fusionForFeature != null) {
                        featuresWithFusionCount++;
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

        LOGGER.info("Analysis performed on {} features in {} VICC entries", totalFeatureCount, resultsPerEntry.size());
        LOGGER.info(" Extracted {} hotspots for {} features", totalHotspotsCount, featuresWithHotspotsCount);
        LOGGER.info(" Extracted {} codons for {} features", totalCodonCount, featuresWithCodonCount);
        LOGGER.info(" Extracted {} exons for {} features", totalExonCount, featuresWithExonCount);
        LOGGER.info(" Extracted {} gene level events", featuresWithGeneLevelEventCount);
        LOGGER.info(" Extracted {} known amps and dels", featuresWithCopyNumberCount);
        LOGGER.info(" Extracted {} known fusions pairs", featuresWithFusionCount);
        LOGGER.info(" Extracted {} signatures", featuresWithSignatureCount);
    }

    public static void writeInterpretationToTsv(@NotNull String featureInterpretationTsv,
            @NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) throws IOException {
        List<String> lines = Lists.newArrayList();
        String header =
                new StringJoiner(FIELD_DELIMITER).add("source").add("gene").add("type").add("feature").add("interpretation").toString();
        lines.add(header);

        for (Map.Entry<ViccEntry, ViccExtractionResult> resultPerEntry : resultsPerEntry.entrySet()) {
            ViccEntry entry = resultPerEntry.getKey();
            ViccExtractionResult result = resultPerEntry.getValue();
            for (Feature feature : entry.features()) {
                StringJoiner interpretation = new StringJoiner(SUB_FIELD_DELIMITER);

                List<VariantHotspot> hotspotsForFeature = result.hotspotsPerFeature().get(feature);
                if (hotspotsForFeature != null) {
                    interpretation.add(hotspotsForFeature.toString());
                }

                List<CodonAnnotation> codonsForFeature = result.codonsPerFeature().get(feature);
                if (codonsForFeature != null) {
                    interpretation.add(codonsForFeature.toString());
                }

                List<ExonAnnotation> exonsForFeature = result.exonsPerFeature().get(feature);
                if (exonsForFeature != null) {
                    interpretation.add(exonsForFeature.toString());
                }

                GeneLevelAnnotation geneLevelEventForFeature = result.geneLevelEventsPerFeature().get(feature);
                if (geneLevelEventForFeature != null) {
                    interpretation.add(geneLevelEventForFeature.toString());
                }

                KnownCopyNumber ampDelForFeature = result.ampsDelsPerFeature().get(feature);
                if (ampDelForFeature != null) {
                    interpretation.add(ampDelForFeature.toString());
                }

                KnownFusionPair fusionForFeature = result.fusionsPerFeature().get(feature);
                if (fusionForFeature != null) {
                    interpretation.add(fusionForFeature.toString());
                }

                SignatureName signatureForFeature = result.signaturesPerFeature().get(feature);
                if (signatureForFeature != null) {
                    interpretation.add(signatureForFeature.toString());
                }

                lines.add(new StringJoiner(FIELD_DELIMITER).add(entry.source().display())
                        .add(feature.geneSymbol())
                        .add(feature.type().toString())
                        .add(feature.name())
                        .add(interpretation.toString())
                        .toString());
            }
        }

        LOGGER.info("Writing {} VICC feature interpretations to {}", lines.size() - 1, featureInterpretationTsv);
        Files.write(new File(featureInterpretationTsv).toPath(), lines);
    }

    public static void writeFeaturesToTsv(@NotNull String featureTsv, @NotNull List<ViccEntry> entries) throws IOException {
        List<String> lines = Lists.newArrayList();
        String header = new StringJoiner(FIELD_DELIMITER).add("source").add("gene").add("transcript").add("type").add("feature").toString();
        lines.add(header);

        Set<FeatureTypeEntry> typeEntries = Sets.newHashSet();

        for (ViccEntry entry : entries) {
            for (Feature feature : entry.features()) {
                typeEntries.add(new FeatureTypeEntry(entry.source().display(),
                        feature.geneSymbol(),
                        entry.transcriptId(),
                        feature.type().toString(),
                        feature.name()));
            }
        }

        for (FeatureTypeEntry feature : typeEntries) {
            lines.add(new StringJoiner(FIELD_DELIMITER).add(feature.source())
                    .add(feature.gene())
                    .add(feature.transcript())
                    .add(feature.type())
                    .add(feature.name())
                    .toString());
        }

        LOGGER.info("Writing {} unique VICC features to {}", lines.size() - 1, featureTsv);
        Files.write(new File(featureTsv).toPath(), lines);
    }

    private static class FeatureTypeEntry {

        @NotNull
        private final String source;
        @Nullable
        private final String gene;
        @Nullable
        private final String transcript;
        @NotNull
        private final String type;
        @NotNull
        private final String name;

        public FeatureTypeEntry(@NotNull final String source, @Nullable final String gene, @Nullable final String transcript,
                @NotNull final String type, @NotNull final String name) {
            this.source = source;
            this.gene = gene;
            this.transcript = transcript;
            this.type = type;
            this.name = name;
        }

        @NotNull
        public String source() {
            return source;
        }

        @Nullable
        public String gene() {
            return gene;
        }

        @Nullable
        public String transcript() {
            return transcript;
        }

        @NotNull
        public String type() {
            return type;
        }

        @NotNull
        public String name() {
            return name;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final FeatureTypeEntry that = (FeatureTypeEntry) o;
            return source.equals(that.source) && Objects.equals(gene, that.gene) && Objects.equals(transcript, that.transcript)
                    && type.equals(that.type) && name.equals(that.name);
        }

        @Override
        public int hashCode() {
            return Objects.hash(source, gene, transcript, type, name);
        }
    }
}
