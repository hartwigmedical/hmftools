package com.hartwig.hmftools.patientreporter.loadStructuralVariants;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class FusionFactory {

    private FusionFactory() {
    }

    private static final String DELIMITER = ",";

    @NotNull
    public static FusionAnalyzer fromFusionFile(@NotNull String fusionFile) throws IOException {
        final List<Fusion> fusions = Lists.newArrayList();

        final List<String> lineFusions = Files.readAllLines(new File(fusionFile).toPath());

        // Skip header line
        for (String line : lineFusions.subList(1, lineFusions.size())) {
            fusions.add(fromFusionLine(line));
        }
        return new FusionAnalyzer(fusions);
    }

    @NotNull
    private static Fusion fromFusionLine(@NotNull String line) {
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
                .proteinsLost(values[12])
                .geneUp(values[13])
                .chrBandUp(values[14])
                .transcriptUp(values[15])
                .strandUp(values[16])
                .regionTypeUp(values[17])
                .codingTypeUp(values[18])
                .exonUp(values[19])
                .phaseUp(values[20])
                .exonMaxUp(values[21])
                .disruptiveUp(values[22])
                .exactBaseUp(values[23])
                .codingBasesUp(values[24])
                .totalCodingUp(values[25])
                .codingStartUp(values[26])
                .codingEndUp(values[27])
                .transStartUp(values[28])
                .transEndUp(values[29])
                .distancePrevUp(values[30])
                .canonicalUp(values[31])
                .biotypeUp(values[32])
                .svIdDown(values[33])
                .chrDown(values[34])
                .posDown(values[35])
                .orientDown(values[36])
                .typeDown(values[37])
                .ploidyDown(values[38])
                .geneDown(values[39])
                .chrBandDown(values[40])
                .transcriptDown(values[41])
                .strandDown(values[42])
                .regionTypeDown(values[43])
                .codingTypeDown(values[44])
                .exonDown(values[45])
                .phaseDown(values[46])
                .exonMaxDown(values[47])
                .disruptiveDown(values[48])
                .exactBaseDown(values[49])
                .codingBasesDown(values[50])
                .totalCodingDown(values[51])
                .codingStartDown(values[52])
                .codingEndDown(values[53])
                .transStartDown(values[54])
                .transEndDown(values[55])
                .distancePrevDown(values[56])
                .canonicalDown(values[57])
                .biotypeDown(values[58])
                .proteinsKept(values.length > 59 ? values[59] : Strings.EMPTY)
                .proteinsLost(values.length > 60 ? values[60] : Strings.EMPTY)
                .build();
    }
}
