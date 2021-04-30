package com.hartwig.hmftools.common.peach;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.jetbrains.annotations.NotNull;

public final class PeachCallsFile {

    private static final String DELIMITER = "\t";

    private PeachCallsFile() {
    }

    @NotNull
    public static List<PeachCalls> read(@NotNull final String filename) throws IOException {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    @NotNull
    private static List<PeachCalls> fromLines(@NotNull final List<String> lines) {
        return lines.stream()
                .skip(1)
                .map(PeachCallsFile::fromString)
                .collect(toList());
    }

    @NotNull
    private static PeachCalls fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);

        final ImmutablePeachCalls.Builder builder = ImmutablePeachCalls.builder()
                .gene(values[0])
                .chromosome(values[1])
                .positionV37(values[2])
                .positionV38(values[3])
                .refV37(values[4])
                .refV38(values[5])
                .allele1(values[6])
                .allele2(values[7])
                .rsid(values[8])
                .variantAnnotationV37(values[9])
                .filterV37(values[10])
                .variantAnnotationV38(values[11])
                .filterV38(values[12])
                .panelVersion(values[13])
                .repoVersion(values[14]);

        return builder.build();
    }
}
