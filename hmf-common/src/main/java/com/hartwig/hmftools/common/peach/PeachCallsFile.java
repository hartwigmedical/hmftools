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
                .positionGRCh37(values[1])
                .refGRCh37(values[2])
                .altGRCh37(values[3])
                .positionGRCh38(values[4])
                .refGRCh38(values[5])
                .altGRCh38(values[6])
                .rsid(values[7])
                .variantAnnotation(values[8])
                .filter(values[9]);

        return builder.build();
    }
}
