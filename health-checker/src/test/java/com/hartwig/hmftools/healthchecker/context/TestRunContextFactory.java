package com.hartwig.hmftools.healthchecker.context;

import java.io.File;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TestRunContextFactory {

    private TestRunContextFactory() {
    }

    @NotNull
    public static RunContext forTest(@NotNull final String runDirectory) {
        return forSomaticTest(runDirectory, Strings.EMPTY, Strings.EMPTY);
    }

    @NotNull
    public static RunContext forSingleSampleTest(@NotNull final String runDirectory, @NotNull final String sample) {
        return new RunContextImpl(runDirectory, removePath(runDirectory), sample, Strings.EMPTY, false, false);
    }

    @NotNull
    public static RunContext forSomaticTest(@NotNull final String runDirectory, @NotNull final String refSample,
            @NotNull final String tumorSample) {
        return new RunContextImpl(runDirectory, removePath(runDirectory), refSample, tumorSample, false, true);
    }

    @NotNull
    private static String removePath(@NotNull final String runDirectory) {
        String folderName = runDirectory;
        if (runDirectory.contains(File.separator)) {
            folderName = runDirectory.substring(runDirectory.lastIndexOf(File.separator) + 1, runDirectory.length());
        }
        return folderName;
    }
}
