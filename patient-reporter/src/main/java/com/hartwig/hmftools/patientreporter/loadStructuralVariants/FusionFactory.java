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
                .ploidyUp(Double.parseDouble(values[12]))
                .proteinsLost(values[13])
                .geneUp(values[14])
                .chrBandUp(values[15])
                .transcriptUp(values[16])
                .strandUp(values[17])
                .regionTypeUp(values[18])
                .codingTypeUp(values[19])
                .exonUp(values[20])
                .phaseUp(values[21])
                .exonMaxUp(values[22])
                .disruptiveUp(values[23])
                .exactBaseUp(values[24])
                .codingBasesUp(values[25])
                .totalCodingUp(values[26])
                .codingStartUp(values[27])
                .codingEndUp(values[28])
                .transStartUp(values[29])
                .transEndUp(values[30])
                .distancePrevUp(values[31])
                .canonicalUp(values[32])
                .biotypeUp(values[33])
                .svIdDown(values[34])
                .chrDown(values[35])
                .posDown(values[36])
                .orientDown(values[37])
                .typeDown(values[38])
                .ploidyDown(Double.parseDouble(values[39]))
                .geneDown(values[40])
                .chrBandDown(values[41])
                .transcriptDown(values[42])
                .strandDown(values[43])
                .regionTypeDown(values[44])
                .codingTypeDown(values[45])
                .exonDown(values[46])
                .phaseDown(values[47])
                .exonMaxDown(values[48])
                .disruptiveDown(values[49])
                .exactBaseDown(values[50])
                .codingBasesDown(values[51])
                .totalCodingDown(values[52])
                .codingStartDown(values[53])
                .codingEndDown(values[54])
                .transStartDown(values[55])
                .transEndDown(values[56])
                .distancePrevDown(values[57])
                .canonicalDown(values[58])
                .biotypeDown(values[59])
                .proteinsKept(values.length > 60 ? values[60] : Strings.EMPTY)
                .proteinsLost(values.length > 61 ? values[61] : Strings.EMPTY)
                .build();
    }
}
