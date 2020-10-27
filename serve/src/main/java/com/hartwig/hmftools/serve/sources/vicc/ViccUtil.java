package com.hartwig.hmftools.serve.sources.vicc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
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
                .add("name")
                .add("type")
                .add("proteinAnnotation")
                .add("transcript")
                .toString();
        lines.add(header);

        for (ViccEntry entry : entries) {
            for (Feature feature : entry.features()) {
                lines.add(new StringJoiner(FIELD_DELIMITER).add(entry.source().display())
                        .add(feature.geneSymbol())
                        .add(feature.name())
                        .add(feature.type().toString())
                        .add(feature.proteinAnnotation())
                        .add(entry.transcriptId())
                        .toString());
            }
        }

        LOGGER.info("Writing {} VICC features to {}", lines.size() - 1, viccFeatureTsv);
        Files.write(new File(viccFeatureTsv).toPath(), lines);
    }
}
