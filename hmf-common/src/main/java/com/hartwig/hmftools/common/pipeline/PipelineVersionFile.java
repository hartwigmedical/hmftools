package com.hartwig.hmftools.common.pipeline;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class PipelineVersionFile {

    private static final Logger LOGGER = LogManager.getLogger(PipelineVersionFile.class);

    private PipelineVersionFile() {
    }

    @Nullable
    public static String majorDotMinorVersion(@NotNull String pipelineVersionFile) throws IOException {
        String version = readPipelineVersion(pipelineVersionFile);
        return version != null ? convertToMajorDotVersion(version) : null;
    }

    @Nullable
    @VisibleForTesting
    static String convertToMajorDotVersion(@NotNull String version) {
        // In case no dot is present there is no major.minor version either.
        if (!version.contains(".")) {
            return null;
        }

        String[] versionSplit = version.split("\\.");
        return versionSplit[0] + "." + versionSplit[1];
    }

    @Nullable
    private static String readPipelineVersion(@NotNull String pipelineVersionFile) throws IOException {
        List<String> lines = Files.readAllLines(new File(pipelineVersionFile).toPath());
        if (lines.isEmpty()) {
            throw new IOException("Pipeline version file seems empty on " + pipelineVersionFile);
        } else {
            if (lines.size() == 1) {
                return lines.get(0);
            } else {
                LOGGER.warn("Too many lines in pipeline version file {}!", pipelineVersionFile);
                return null;
            }
        }
    }
}
