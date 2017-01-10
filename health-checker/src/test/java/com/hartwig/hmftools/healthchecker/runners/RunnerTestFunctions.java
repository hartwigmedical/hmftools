package com.hartwig.hmftools.healthchecker.runners;

import java.io.File;

import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;

final class RunnerTestFunctions {

    private RunnerTestFunctions() {
    }

    @NotNull
    static String getRunnerResourcePath(@NotNull String runner) {
        return Resources.getResource("Runners" + File.separator + runner).getPath();
    }
}
