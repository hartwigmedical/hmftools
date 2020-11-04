package com.hartwig.hmftools.serve.sources.vicc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusionFile;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGeneFile;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspotFile;
import com.hartwig.hmftools.serve.actionability.range.ActionableRangeFile;
import com.hartwig.hmftools.serve.actionability.signature.ActionableSignatureFile;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ViccUtil {

    private static final Logger LOGGER = LogManager.getLogger(ViccUtil.class);

    private static final String FIELD_DELIMITER = "\t";

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

    public static void writeFeatures(@NotNull String viccFeatureTsv, @NotNull List<ViccEntry> entries) throws IOException {
        List<String> lines = Lists.newArrayList();
        String header = new StringJoiner(FIELD_DELIMITER).add("source")
                .add("gene")
                .add("transcript")
                .add("type")
                .add("name")
                .toString();
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

        LOGGER.info("Writing {} unique VICC features to {}", lines.size() - 1, viccFeatureTsv);
        Files.write(new File(viccFeatureTsv).toPath(), lines);
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
