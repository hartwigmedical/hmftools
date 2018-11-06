package com.hartwig.hmftools.patientreporter.chord;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ImmutableChordAnalysis;

import org.jetbrains.annotations.NotNull;

public final class ChordFile {
    private static final String DELIMITER = "\t";
    private static final String CHORD_FILE_EXTENSION = "_chord_prediction.txt";

    private ChordFile() {
    }

    @NotNull
    public static ChordAnalysis loadChordFile(@NotNull String chordDirectory, @NotNull String sample) throws IOException {
        return fromFile(chordDirectory + File.separator + sample + CHORD_FILE_EXTENSION);
    }

    @VisibleForTesting
    @NotNull
    static ChordAnalysis fromFile(@NotNull String filePath) throws IOException {
        List<String> allLines = Files.readAllLines(new File(filePath).toPath());
        assert allLines.size() == 2;

        return fromLine(allLines.get(1));
    }

    @NotNull
    private static ChordAnalysis fromLine(@NotNull String line) {
        String[] values = line.split(DELIMITER);

        return ImmutableChordAnalysis.builder()
                .noneValue(Double.valueOf(values[2]))
                .BRCA1Value(Double.valueOf(values[1]))
                .BRCA2Value(Double.valueOf(values[3]))
                .hrdValue(Double.valueOf(values[4]))
                .predictedResponseValue(Double.valueOf(values[5]))
                .build();
    }
}
