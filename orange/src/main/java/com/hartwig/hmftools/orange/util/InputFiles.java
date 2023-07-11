package com.hartwig.hmftools.orange.util;

import static java.lang.String.format;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class InputFiles {

    private InputFiles() {
    }

    public static String directoryExists(@NotNull String directory) {
        if (!Files.isDirectory(Path.of(directory))) {
            throw new IllegalArgumentException(format("Path [%s] is not a directory", directory));
        }
        return fileExists(directory);
    }

    public static String fileExists(@NotNull String fullPath) {
        boolean exists = Files.exists(new File(fullPath).toPath());
        if (!exists) {
            throw new IllegalArgumentException(format("Path [%s] did not exist and is mandatory.", fullPath));
        }
        return fullPath;
    }
}
