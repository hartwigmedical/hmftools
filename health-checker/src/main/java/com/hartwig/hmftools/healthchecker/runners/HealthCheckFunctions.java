package com.hartwig.hmftools.healthchecker.runners;

import static java.util.stream.Collectors.toList;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.healthchecker.exception.EmptyFileException;
import com.hartwig.hmftools.healthchecker.exception.HealthChecksException;
import com.hartwig.hmftools.healthchecker.io.reader.ZipFilesReader;

import org.jetbrains.annotations.NotNull;

final class HealthCheckFunctions {

    private static final String FASTQC_DATA_FILE_NAME = "fastqc_data.txt";
    private static final String TOTAL_SEQUENCES_PATTERN = "Total Sequences";
    private static final String COLUMN_SEPARATOR = "\t";

    private HealthCheckFunctions() {
    }

    static long sumOfTotalSequencesFromFastQC(@NotNull final String basePath)
            throws IOException, HealthChecksException {
        final List<String> allLines = new ZipFilesReader().readFieldFromZipFiles(basePath, FASTQC_DATA_FILE_NAME,
                TOTAL_SEQUENCES_PATTERN);

        final List<String> allValues = allLines.stream().map(line -> {
            String totalSequences = null;
            if (line != null) {
                final String[] values = line.split(COLUMN_SEPARATOR);
                totalSequences = values[1];
            }
            return totalSequences;
        }).filter(lines -> lines != null).collect(toList());

        long totalSequences = !allValues.isEmpty() ? allValues.stream().mapToLong(Long::parseLong).sum() : 0L;
        if (totalSequences == 0) {
            throw new EmptyFileException(FASTQC_DATA_FILE_NAME, basePath);
        }
        return totalSequences;
    }
}
