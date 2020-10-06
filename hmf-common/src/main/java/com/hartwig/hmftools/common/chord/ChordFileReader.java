package com.hartwig.hmftools.common.chord;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Optional;
import java.util.stream.IntStream;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.utils.io.exception.EmptyFileException;
import com.hartwig.hmftools.common.utils.io.exception.MalformedFileException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ChordFileReader {

    private static final String VALUE_SEPARATOR = "\t";
    @VisibleForTesting
    static final String V1_NA = "UNKNOWN";

    private static final int V1_BRCA1_COLUMN = 2;
    private static final int V1_BRCA2_COLUMN = 3;
    private static final int V1_HRD_COLUMN = 4;

    private static final int V2_BRCA1_COLUMN = 1;
    private static final int V2_BRCA2_COLUMN = 2;
    private static final int V2_HRD_COLUMN = 3;
    private static final int V2_HR_STATUS_COLUMN = 4;
    private static final int V2_HRD_TYPE_COLUMN = 5;
    private static final int V2_REMARKS_HR_STATUS_COLUMN = 6;
    private static final int V2_REMARKS_HRD_TYPE_COLUMN = 7;

    private static final Logger LOGGER = LogManager.getLogger(ChordFileReader.class);


    private ChordFileReader() {
    }

    @NotNull
    public static ChordAnalysis read(@NotNull String filePath) throws IOException {
        ImmutableChordAnalysis.Builder builder = ImmutableChordAnalysis.builder();
        String[] values = findValuesLine(filePath).split(VALUE_SEPARATOR);

        if (isV1(values)) {
            builder.BRCA1Value(Double.parseDouble(values[V1_BRCA1_COLUMN]));
            builder.BRCA2Value(Double.parseDouble(values[V1_BRCA2_COLUMN]));
            builder.hrdValue(Double.parseDouble(values[V1_HRD_COLUMN]));
            builder.hrStatus(extractHrStatus(V1_NA));
            builder.hrdType(V1_NA);
            builder.remarksHrStatus(V1_NA);
            builder.remarksHrdType(V1_NA);
        } else {
            String remarksHrStatus = values.length > V2_REMARKS_HR_STATUS_COLUMN ? values[V2_REMARKS_HR_STATUS_COLUMN] : Strings.EMPTY;
            String remarksHrdType = values.length > V2_REMARKS_HRD_TYPE_COLUMN ? values[V2_REMARKS_HRD_TYPE_COLUMN] : Strings.EMPTY;

            builder.BRCA1Value(Double.parseDouble(values[V2_BRCA1_COLUMN]));
            builder.BRCA2Value(Double.parseDouble(values[V2_BRCA2_COLUMN]));
            builder.hrdValue(Double.parseDouble(values[V2_HRD_COLUMN]));
            builder.hrStatus(extractHrStatus(values[V2_HR_STATUS_COLUMN]));
            builder.hrdType(values[V2_HRD_TYPE_COLUMN]);
            builder.remarksHrStatus(remarksHrStatus);
            builder.remarksHrdType(remarksHrdType);
        }

        return builder.build();
    }

    @NotNull
    public static ChordStatus extractHrStatus(@NotNull String hrStatus) {
        switch (hrStatus) {
            case "cannot_be_determined": return ChordStatus.CANNOT_BE_DETERMINED;
            case "HR_proficient": return ChordStatus.HR_PROFICIENT;
            case "HR_deficient": return ChordStatus.HR_DEFICIENT;
        }
        LOGGER.warn("Unknown HR status of CHORD");
        return ChordStatus.UNKNOWN;
    }

    @NotNull
    private static String findValuesLine(@NotNull String filename) throws IOException {
        List<String> lines = Files.readAllLines(new File(filename).toPath());
        if (lines.isEmpty()) {
            throw new EmptyFileException(filename);
        }
        int index = findHeaderLineIndex(lines);
        if (index >= lines.size()) {
            throw new MalformedFileException(String.format("No value line found after header line in CHORD file %s.", filename));
        }
        return lines.get(index + 1);
    }

    private static int findHeaderLineIndex(@NotNull List<String> lines) throws MalformedFileException {
        Optional<Integer> lineNumbers =
                IntStream.range(0, lines.size()).filter(index -> lines.get(index).contains("hrd")).boxed().findFirst();
        if (!lineNumbers.isPresent()) {
            throw new MalformedFileException(String.format("Could not find header line in CHORD file with %s lines.", lines.size()));
        }
        return lineNumbers.get();
    }

    private static boolean isV1(@NotNull String[] values) {
        assert V1_HRD_COLUMN == V2_HR_STATUS_COLUMN && V1_HRD_COLUMN == 4;
        return isDouble(values[4]);
    }

    private static boolean isDouble(@NotNull String value) {
        try {
            Double.parseDouble(value);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }
}
