package com.hartwig.hmftools.orange.report;

import java.io.File;
import java.nio.file.Files;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PlotPathResolver {

    @Nullable
    private final String outputDir;

    PlotPathResolver(@Nullable final String outputDir) {
        this.outputDir = outputDir;
    }

    @NotNull
    public String resolve(@NotNull String path) {
        if (outputDir == null || Files.exists(new File(path).toPath())) {
            return path;
        }

        // Assume if path in absolute does not exist, the path is relative to the outputDir.
        return outputDir + File.separator + path;
    }
}
