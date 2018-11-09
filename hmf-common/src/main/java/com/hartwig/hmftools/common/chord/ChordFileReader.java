package com.hartwig.hmftools.common.chord;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Optional;
import java.util.stream.IntStream;

import com.hartwig.hmftools.common.io.exception.EmptyFileException;
import com.hartwig.hmftools.common.io.exception.MalformedFileException;
import com.hartwig.hmftools.common.io.path.PathPrefixSuffixFinder;

import org.jetbrains.annotations.NotNull;

public final class ChordFileReader {

    // LISC: chord files stores in {run}/chord/{sample}_chord_prediction.txt
    private static final String CHORD_BASE_DIRECTORY = "chord";
    private static final String CHORD_EXTENSION = "_chord_prediction.txt";
    private static final String VALUE_SEPARATOR = "\t";

    private static final int BRCA1_COLUMN = 1;
    private static final int NONE_COLUMN = 2;
    private static final int BRCA2_COLUMN = 3;
    private static final int HRD_COLUMN = 4;
    private static final int PREDICTED_RESPONSE_COLUMN = 5;

    private ChordFileReader() {
    }

    @NotNull
    public static String generateFilename(@NotNull final String runDir, @NotNull final String sample) throws FileNotFoundException {
        String path = runDir + File.separator + CHORD_BASE_DIRECTORY;
        return PathPrefixSuffixFinder.build().findPath(path, sample, CHORD_EXTENSION).toString();
    }

    @NotNull
    public static ChordAnalysis read(@NotNull String filePath) throws IOException {
        ImmutableChordAnalysis.Builder builder = ImmutableChordAnalysis.builder();

        builder = appendChordValues(builder, filePath);

        return builder.build();
    }

    @NotNull
    private static ImmutableChordAnalysis.Builder appendChordValues(@NotNull ImmutableChordAnalysis.Builder builder,
            @NotNull String filePath) throws IOException {
        String[] values = findValuesLine(filePath).split(VALUE_SEPARATOR);
        builder.BRCA1Value(Double.valueOf(values[BRCA1_COLUMN]));
        builder.noneValue(Double.valueOf(values[NONE_COLUMN]));
        builder.BRCA2Value(Double.valueOf(values[BRCA2_COLUMN]));
        builder.hrdValue(Double.valueOf(values[HRD_COLUMN]));
        builder.predictedResponseValue(Double.valueOf(values[PREDICTED_RESPONSE_COLUMN]));
        return builder;
    }

    @NotNull
    private static String findValuesLine(@NotNull String filename) throws IOException {
        List<String> lines = Files.readAllLines(new File(filename).toPath());
        if (lines.isEmpty()) {
            throw new EmptyFileException(filename);
        }
        final int index = findHeaderLineIndex(lines);
        if (index >= lines.size()) {
            throw new MalformedFileException(String.format("No value line found after header line in CHORD file %s.", filename));
        }
        return lines.get(index + 1);
    }

    private static int findHeaderLineIndex(@NotNull final List<String> lines) throws MalformedFileException {
        final Optional<Integer> lineNumbers =
                IntStream.range(0, lines.size()).filter(index -> lines.get(index).contains("hrd")).boxed().findFirst();
        if (!lineNumbers.isPresent()) {
            throw new MalformedFileException(String.format("Could not find header line in CHORD file with %s lines.", lines.size()));
        }
        return lineNumbers.get();
    }

}
