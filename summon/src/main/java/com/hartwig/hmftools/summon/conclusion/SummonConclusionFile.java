package com.hartwig.hmftools.summon.conclusion;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class SummonConclusionFile {

    private static final String EXTENSION = ".summon.tsv";

    private SummonConclusionFile() {
    }

    @NotNull
    public static String generateFilename(@NotNull String basePath, @NotNull String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static void write(@NotNull String file, @NotNull ActionabilityConclusion actionabilityConclusion) throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(toLine(actionabilityConclusion.conclusion()));
        Files.write(new File(file).toPath(), lines);
    }

    @NotNull
    private static String toLine(@NotNull String conclusion) {
        return conclusion;
    }
}