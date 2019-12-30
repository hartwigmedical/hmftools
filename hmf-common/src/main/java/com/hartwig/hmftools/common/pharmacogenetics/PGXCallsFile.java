package com.hartwig.hmftools.common.pharmacogenetics;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.common.utils.io.exception.MalformedFileException;

import org.jetbrains.annotations.NotNull;

public final class PGXCallsFile {

    private PGXCallsFile() {
    }


    @NotNull
    public static PGXCalls read(@NotNull final String filename) throws IOException {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    @NotNull
    private static PGXCalls fromLines(@NotNull final List<String> lines) throws MalformedFileException {
        try {
            return ImmutablePGXCalls.builder()
                    .gene(lines.get(0))
                    .positionGRCh37(lines.get(1))
                    .refGRCh37(lines.get(2))
                    .altGRCh37(lines.get(3))
                    .positionGRCh38(lines.get(4))
                    .refGRCh38(lines.get(5))
                    .altGRCh38(lines.get(6))
                    .rsid(lines.get(7))
                    .variantAnnotation(lines.get(8))
                    .filter(lines.get(9))
                    .build();
        } catch (Exception e) {
            throw new MalformedFileException("Unable to parse purple qc file.");
        }

    }
}
