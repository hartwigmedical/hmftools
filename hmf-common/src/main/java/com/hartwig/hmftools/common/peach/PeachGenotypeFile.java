package com.hartwig.hmftools.common.peach;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.jetbrains.annotations.NotNull;

public final class PeachGenotypeFile {

    private static final String DELIMITER = "\t";

    private PeachGenotypeFile() {
    }

    @NotNull
    public static List<PeachGenotype> read(@NotNull final String filename) throws IOException {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    @NotNull
    private static List<PeachGenotype> fromLines(@NotNull final List<String> lines) {
        return lines.stream()
                .skip(1)
                .map(PeachGenotypeFile::fromString)
                .collect(toList());
    }

    @NotNull
    private static PeachGenotype fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);

        final ImmutablePeachGenotype.Builder builder = ImmutablePeachGenotype.builder()
                .gene(values[0])
                .haplotype(values[1])
                .function(values[2])
                .linkedDrugs(values[3])
                .urlPrescriptionInfo(values[4])
                .panelVersion(values[5])
                .repoVersion(values[6]);

        return builder.build();
    }
}
