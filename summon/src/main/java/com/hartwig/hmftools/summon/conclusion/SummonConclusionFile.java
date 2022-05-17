package com.hartwig.hmftools.summon.conclusion;

import java.io.File;

import org.jetbrains.annotations.NotNull;

public class SummonConclusionFile {

    private static final String EXTENSION = ".summon.tsv";

    private SummonConclusionFile(){}

    @NotNull
    public static String generateFilename(@NotNull String basePath, @NotNull String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }
}