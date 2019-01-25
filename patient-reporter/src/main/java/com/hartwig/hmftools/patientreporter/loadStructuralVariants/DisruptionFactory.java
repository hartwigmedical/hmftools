package com.hartwig.hmftools.patientreporter.loadStructuralVariants;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class DisruptionFactory {

    private DisruptionFactory() {
    }

    private static final String DELIMITER = ",";

    @NotNull
    public static DisruptionAnalyzer fromDisruptionFile(@NotNull String disruptionFile) throws IOException {
        final List<Disruption> disruptions = Lists.newArrayList();

        final List<String> lineDisruptions = Files.readAllLines(new File(disruptionFile).toPath());

        // Skip header line
        for (String line : lineDisruptions.subList(1, lineDisruptions.size())) {
            disruptions.add(fromDisruptionLine(line));
        }
        return new DisruptionAnalyzer(disruptions);
    }

    @NotNull
    private static Disruption fromDisruptionLine(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableDisruption.builder()
                .reportable(Boolean.valueOf(values[1]))
                .svId(values[2])
                .chromosome(values[3])
                .position(values[4])
                .orientation(values[5])
                .type(values[6])
                .ploidy(Double.valueOf(values[7]))
                .gene(values[8])
                .chrBand(values[9])
                .transcript(values[10])
                .strand(values[11])
                .regionType(values[12])
                .codingType(values[13])
                .canonical(values[14])
                .biotype(values[15])
                .exonUp(values[16])
                .exonDown(values[17])
                .isDisruptive(Boolean.valueOf(values[18]))
                .build();
    }
}
