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
                .gene(values[7])
                .transcript(values[8])
                .strand(values[9])
                .regionType(values[10])
                .codingType(values[11])
                .biotype(values[12])
                .exon(values[13])
                .isDisruptive(Boolean.valueOf(values[14]))
                .build();
    }
}
