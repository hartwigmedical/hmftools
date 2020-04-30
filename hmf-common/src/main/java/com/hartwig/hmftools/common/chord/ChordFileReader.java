package com.hartwig.hmftools.common.chord;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Optional;
import java.util.stream.IntStream;

import com.hartwig.hmftools.common.utils.io.exception.EmptyFileException;
import com.hartwig.hmftools.common.utils.io.exception.MalformedFileException;

import org.jetbrains.annotations.NotNull;

public final class ChordFileReader {

    private static final String VALUE_SEPARATOR = "\t";

    private static final int NONE_COLUMN = 1;
    private static final int BRCA1_COLUMN = 2;
    private static final int BRCA2_COLUMN = 3;
    private static final int HRD_COLUMN = 4;
    private static final int PREDICTED_RESPONSE_COLUMN = 5;

    private ChordFileReader() {
    }

    @NotNull
    public static ChordAnalysis read(@NotNull String filePath) throws IOException {
        ImmutableChordAnalysis.Builder builder = ImmutableChordAnalysis.builder();
        String[] values = findValuesLine(filePath).split(VALUE_SEPARATOR);

        builder.noneValue(Double.parseDouble(values[NONE_COLUMN]));
        builder.BRCA1Value(Double.parseDouble(values[BRCA1_COLUMN]));
        builder.BRCA2Value(Double.parseDouble(values[BRCA2_COLUMN]));
        builder.hrdValue(Double.parseDouble(values[HRD_COLUMN]));
        builder.predictedResponseValue(Boolean.parseBoolean(values[PREDICTED_RESPONSE_COLUMN]));

        return builder.build();
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
}
