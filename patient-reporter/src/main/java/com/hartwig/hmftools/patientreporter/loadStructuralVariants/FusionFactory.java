package com.hartwig.hmftools.patientreporter.loadStructuralVariants;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class FusionFactory {

    private FusionFactory() {
    }

    private static final String DELIMITER = ",";

    @NotNull
    public static FusionAnalyzer readingFusion(@NotNull String fusionFile)
            throws IOException {

        final List<Fusion> fusions = new ArrayList<>();

        final List<String> lineFusions = Files.readAllLines(new File(fusionFile).toPath());

        // Skip header line
        for (String line : lineFusions.subList(1, lineFusions.size())) {
            fusions.add(fromLineVariants(line));
        }
        return new FusionAnalyzer(fusions);
    }

    @NotNull
    private static Fusion fromLineVariants(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        // ProteinsKept and ProteinsLost are not mandatory
        return ImmutableFusion.builder()
                .reportable(Boolean.valueOf(values[1]))
                .knownType(values[2])
                .primarySource(values[3])
                .clusterId(values[4])
                .clusterCount(values[5])
                .resolvedType(values[6])
                .svIdUp(values[7])
                .chrUp(values[8])
                .posUp(values[9])
                .orientUp(values[10])
                .typeUp(values[11])
                .geneUp(values[12])
                .transcriptUp(values[13])
                .strandUp(values[14])
                .regionTypeUp(values[15])
                .codingTypeUp(values[16])
                .exonUp(values[17])
                .phaseUp(values[18])
                .exonMaxUp(values[19])
                .disruptiveUp(values[20])
                .exactBaseUp(values[21])
                .codingBasesUp(values[22])
                .totalCodingUp(values[23])
                .codingStartUp(values[24])
                .codingEndUp(values[25])
                .transStartUp(values[26])
                .transEndUp(values[27])
                .distancePrevUp(values[28])
                .biotypeUp(values[29])
                .svIdDown(values[30])
                .chrDown(values[31])
                .posDown(values[32])
                .orientDown(values[33])
                .typeDown(values[34])
                .geneDown(values[35])
                .transcriptDown(values[36])
                .strandDown(values[37])
                .regionTypeDown(values[38])
                .codingTypeDown(values[39])
                .exonDown(values[40])
                .phaseDown(values[41])
                .exonMaxDown(values[42])
                .disruptiveDown(values[43])
                .exactBaseDown(values[44])
                .codingBasesDown(values[45])
                .totalCodingDown(values[46])
                .codingStartDown(values[47])
                .codingEndDown(values[48])
                .transStartDown(values[49])
                .transEndDown(values[50])
                .distancePrevDown(values[51])
                .biotypeDown(values[52])
                .proteinsKept(values.length > 53 ? values[53] : Strings.EMPTY )
                .proteinsLost(values.length > 54 ? values[54] : Strings.EMPTY)
                .build();
    }
}
