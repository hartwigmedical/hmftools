package com.hartwig.hmftools.common.chord;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Optional;
import java.util.stream.IntStream;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ChordFileReader {

    private static final Logger LOGGER = LogManager.getLogger(ChordFileReader.class);

    private static final String VALUE_SEPARATOR = "\t";

    private static final int BRCA1_COLUMN = 1;
    private static final int BRCA2_COLUMN = 2;
    private static final int HRD_COLUMN = 3;
    private static final int HR_STATUS_COLUMN = 4;
    private static final int HRD_TYPE_COLUMN = 5;
    private static final int REMARKS_HR_STATUS_COLUMN = 6;
    private static final int REMARKS_HRD_TYPE_COLUMN = 7;

    private ChordFileReader() {
    }

    @NotNull
    public static ChordAnalysis read(@NotNull String filePath) throws IOException {
        String[] values = findValuesLine(filePath).split(VALUE_SEPARATOR);

        String remarksHrStatus = values.length > REMARKS_HR_STATUS_COLUMN ? values[REMARKS_HR_STATUS_COLUMN] : Strings.EMPTY;
        String remarksHrdType = values.length > REMARKS_HRD_TYPE_COLUMN ? values[REMARKS_HRD_TYPE_COLUMN] : Strings.EMPTY;

        return ImmutableChordAnalysis.builder()
                .BRCA1Value(Double.parseDouble(values[BRCA1_COLUMN]))
                .BRCA2Value(Double.parseDouble(values[BRCA2_COLUMN]))
                .hrdValue(Double.parseDouble(values[HRD_COLUMN]))
                .hrStatus(extractHrStatus(values[HR_STATUS_COLUMN]))
                .hrdType(values[HRD_TYPE_COLUMN])
                .remarksHrStatus(remarksHrStatus)
                .remarksHrdType(remarksHrdType)
                .build();
    }

    @VisibleForTesting
    @NotNull
    static ChordStatus extractHrStatus(@NotNull String hrStatus) {
        switch (hrStatus) {
            case "cannot_be_determined":
                return ChordStatus.CANNOT_BE_DETERMINED;
            case "HR_proficient":
                return ChordStatus.HR_PROFICIENT;
            case "HR_deficient":
                return ChordStatus.HR_DEFICIENT;
        }
        LOGGER.warn("Unknown CHORD HR status: '{}'", hrStatus);
        return ChordStatus.UNKNOWN;
    }

    @NotNull
    private static String findValuesLine(@NotNull String filename) throws IOException {
        List<String> lines = Files.readAllLines(new File(filename).toPath());
        if (lines.isEmpty()) {
            throw new IOException("CHORD file seems empty on " + filename);
        }
        int index = findHeaderLineIndex(lines);
        if (index >= lines.size()) {
            throw new IOException(String.format("No value line found after header line in CHORD file %s.", filename));
        }
        return lines.get(index + 1);
    }

    private static int findHeaderLineIndex(@NotNull List<String> lines) throws IOException {
        Optional<Integer> lineNumbers =
                IntStream.range(0, lines.size()).filter(index -> lines.get(index).contains("hrd")).boxed().findFirst();
        if (!lineNumbers.isPresent()) {
            throw new IOException(String.format("Could not find header line in CHORD file with %s lines.", lines.size()));
        }
        return lineNumbers.get();
    }
}
