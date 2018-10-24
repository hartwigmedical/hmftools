package com.hartwig.hmftools.patientreporter.chordclassifier;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ChordFile {
    private static final Logger LOGGER = LogManager.getLogger(ChordFile.class);
    private static final String DELIMITER = ",";
    private static final String CHORD_FILE_EXTENSION = "_germline_variants.csv";
    private static final String CHORD_HAS_RUN_FILE = "bachelor_found_no_variants";

    private ChordFile() {
    }

    public static boolean hasChordRun(@NotNull String hrdDirectory, @NotNull String sample) {
        if (findChordFilePath(hrdDirectory, sample) == null) {
            String hrdHasRunPath = hrdDirectory + File.separator + CHORD_HAS_RUN_FILE;
            return Files.exists(new File(hrdHasRunPath).toPath());
        }

        return true;
    }

    @Nullable
    public static String findChordFilePath(@NotNull String hrdDirectory, @NotNull String sample) {
        String path = hrdDirectory + File.separator + sample + CHORD_FILE_EXTENSION;
        return Files.exists(new File(path).toPath()) ? path : null;
    }



    @NotNull
    public static List<ChordAnalysis> loadChordFile(@NotNull String filePath) throws IOException {
        Path path = new File(filePath).toPath();
        assert Files.exists(path);
        return fromLines(Files.readAllLines(path));
    }

    @NotNull
    private static List<ChordAnalysis> fromLines(@NotNull List<String> lines) {
        List<ChordAnalysis> chordValues = Lists.newArrayList();
        // KODU: Skip header line
        for (String line : lines.subList(1, lines.size())) {
            chordValues.add(fromString(line));
        }
        return chordValues;
    }

    @NotNull
    private static ChordAnalysis fromString(@NotNull String line) {
        String[] values = line.split(DELIMITER);

        String program = values[1];
        if (!program.equalsIgnoreCase("lynparza")) {
            LOGGER.warn("Unexpected file for chord detecting found: " + program);
        }

        return ImmutableChordAnalysis.builder()
                .BRCA1Value(Double.valueOf(values[0]))
                .noneValue(Double.valueOf(values[1]))
                .BRCA2Value(Double.valueOf(values[2]))
                .hrdValue(Double.valueOf(values[3]))
                .predictedResponseValue(Double.valueOf(values[4]))
                .build();
    }

}
