package com.hartwig.hmftools.protect.conclusion;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class TumorLocationConclusionFile {

    private static final String DELIMITER = "\t";

    private TumorLocationConclusionFile(){

    }

    @NotNull
    public static List<TumorLocationConclusion> readTumorLocationConclusion(@NotNull String readTumorLocationConclusion) throws IOException {
        final List<TumorLocationConclusion> tumorLocationConclusion = Lists.newArrayList();

        final List<String> lineTumorLocationConclusion = Files.readAllLines(new File(readTumorLocationConclusion).toPath());
        for (String tumorConclusion : lineTumorLocationConclusion.subList(1, lineTumorLocationConclusion.size())) {
            tumorLocationConclusion.add(lineTumorLocationConclusion(tumorConclusion));
        }
        return tumorLocationConclusion;
    }

    @NotNull
    private static TumorLocationConclusion lineTumorLocationConclusion(@NotNull String line) {
        final String[] values = line.split(DELIMITER);

       return ImmutableTumorLocationConclusion.builder()
               .primaryTumorLocation(values[0])
               .cancerSubType(values[1])
               .tumorLocationConclusion(values[2])
               .build();
    }
}
