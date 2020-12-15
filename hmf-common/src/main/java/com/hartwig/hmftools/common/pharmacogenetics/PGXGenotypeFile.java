package com.hartwig.hmftools.common.pharmacogenetics;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.jetbrains.annotations.NotNull;

public final class PGXGenotypeFile {

    private static final String DELIMITER = "\t";

    private PGXGenotypeFile() {
    }

    @NotNull
    public static List<PGXGenotype> read(@NotNull final String filename) throws IOException {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    @NotNull
    private static List<PGXGenotype> fromLines(@NotNull final List<String> lines) {
        return lines.stream()
                .skip(1)
                .map(PGXGenotypeFile::fromString)
                .collect(toList());
    }

    @NotNull
    private static PGXGenotype fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);

        final ImmutablePGXGenotype.Builder builder = ImmutablePGXGenotype.builder()
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
