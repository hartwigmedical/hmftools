package com.hartwig.hmftools.common.pharmacogenetics;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.common.utils.io.exception.MalformedFileException;

import org.jetbrains.annotations.NotNull;

public final class PGXGenotypeFile {

    private PGXGenotypeFile() {
    }


    @NotNull
    public static PGXGenotype read(@NotNull final String filename) throws IOException {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    @NotNull
    private static PGXGenotype fromLines(@NotNull final List<String> lines) throws MalformedFileException {
        try {
            return ImmutablePGXGenotype.builder()
                    .gene(lines.get(0))
                    .haplotype(lines.get(1))
                    .function(lines.get(2))
                    .linkedDrugs(lines.get(3))
                    .urlPrescriptionInfo(lines.get(4))
                    .panelVersion(lines.get(5))
                    .repoVersion(lines.get(6))
                    .build();
        } catch (Exception e) {
            throw new MalformedFileException("Unable to parse purple qc file.");
        }
    }

}
