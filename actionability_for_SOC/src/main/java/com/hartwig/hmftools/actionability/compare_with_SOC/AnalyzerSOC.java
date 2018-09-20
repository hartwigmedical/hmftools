package com.hartwig.hmftools.actionability.compare_with_SOC;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.jetbrains.annotations.NotNull;

public class AnalyzerSOC {
    private static final org.apache.logging.log4j.Logger LOGGER = LogManager.getLogger(AnalyzerSOC.class);
    private static final String DELIMITER = "\t";

    @NotNull
    private final List<ReadingRegionsBedFile> dataFile;

    public AnalyzerSOC(@NotNull final List<ReadingRegionsBedFile> dataFile) {
        this.dataFile = dataFile;
    }

    @NotNull
    public static AnalyzerSOC loadFileBedFile(String file) throws IOException {
        final List<ReadingRegionsBedFile> dataFile = new ArrayList<>();
        final List<String> line = Files.readAllLines(new File(file).toPath());

        for (int i = 1; i < line.size(); i++) {
            fromLine(line.get(i));
            dataFile.add(fromLine(line.get(i)));
        }
        LOGGER.info(dataFile);
        return new AnalyzerSOC(dataFile);
    }

    @NotNull
    private static ReadingRegionsBedFile fromLine(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableReadingRegionsBedFile.builder()
                .chromosome(values[0])
                .startPosition(values[1])
                .endposition(values[2])
                .gene(values[3])
                .strand(values[5])
                .build();
    }
}
