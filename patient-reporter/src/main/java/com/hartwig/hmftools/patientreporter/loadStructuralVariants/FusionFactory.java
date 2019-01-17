package com.hartwig.hmftools.patientreporter.loadStructuralVariants;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.jetbrains.annotations.NotNull;

public class FusionFactory {

    private FusionFactory() {
    }

    private static final String DELIMITER = ",";


    @NotNull
    public static FusionAnalyzer readingFusion(@NotNull String fusionFile) throws IOException {
        final List<FusionReaderFile> fusions = new ArrayList<>();

        final List<String> lineFusions = Files.readAllLines(new File(fusionFile).toPath());

        // Skip header line
        for (String linefusion : lineFusions.subList(1, lineFusions.size())) {
            fusions.add(fromLineVariants(linefusion));
        }
        return new FusionAnalyzer (fusions);
    }

    @NotNull
    private static FusionReaderFile fromLineVariants(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableFusionReaderFile.builder()
                .reportable(Boolean.valueOf(values[0]))
                .knownType(values[1])
                .primarySource(values[2])
                .clusterId(values[3])
                .clusterCount(values[4])
                .resolvedType(values[5])
                .svIdUp(values[6])
                .chrUp(values[7])
                .posUp(values[8])
                .orientUp(values[9])
                .typeUp(values[10])
                .geneUp(values[11])
                .transcriptUp(values[12])
                .strandUp(values[13])
                .regionTypeUp(values[14])
                .codingTypeUp(values[15])
                .exonUp(values[16])
                .phaseUp(values[17])
                .exonMaxUp(values[18])
                .disruptiveUp(values[19])
                .exactBaseUp(values[20])
                .codingBasesUp(values[21])
                .totalCodingUp(values[22])
                .codingStartUp(values[23])
                .codingEndUp(values[24])
                .transStartUp(values[25])
                .transEndUp(values[26])
                .distancePrevUp(values[27])
                .biotypeUp(values[28])
                .svIdDown(values[29])
                .chrDown(values[30])
                .posDown(values[31])
                .orientDown(values[32])
                .typeDown(values[33])
                .geneDown(values[34])
                .transcriptDown(values[35])
                .strandDown(values[36])
                .regionTypeDown(values[37])
                .codingTypeDown(values[38])
                .exonDown(values[39])
                .phaseDown(values[40])
                .exonMaxDown(values[41])
                .disruptiveDown(values[42])
                .exactBaseDown(values[43])
                .codingBasesDown(values[44])
                .totalCodingDown(values[45])
                .codingStartDown(values[46])
                .codingEndDown(values[47])
                .transStartDown(values[48])
                .transEndDown(values[49])
                .distancePrevDown(values[50])
                .biotypeDown(values[51])
                .proteinsKept(values[52])
                .proteinsLost(values[53])
                .build();
    }
}
