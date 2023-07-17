package com.hartwig.hmftools.orange.util;

import static java.lang.String.format;

import java.io.File;
import java.nio.file.Files;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class Config {

    private Config() {
    }

    @Nullable
    public static String optionalFileIfExists(@NotNull String fullPath) {
        boolean exists = Files.exists(new File(fullPath).toPath());
        return exists ? fullPath : null;
    }

    @NotNull
    public static String fileIfExists(@NotNull String fullPath) {
        boolean exists = Files.exists(new File(fullPath).toPath());
        if (!exists) {
            throw new IllegalArgumentException(format("Path [%s] did not exist and is mandatory.", fullPath));
        }
        return fullPath;
    }
}
