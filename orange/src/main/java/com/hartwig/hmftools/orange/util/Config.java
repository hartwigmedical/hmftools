package com.hartwig.hmftools.orange.util;

import static java.lang.String.format;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public final class Config {

    private Config() {
    }

    @NotNull
    public static String nonOptionalValue(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);
        if (value == null) {
            throw new IllegalArgumentException("Parameter must be provided: " + param);
        }

        return value;
    }

    @NotNull
    public static String outputDir(@NotNull CommandLine cmd, @NotNull String param) throws ParseException, IOException {
        String value = nonOptionalValue(cmd, param);
        File outputDir = new File(value);
        if (!outputDir.exists() && !outputDir.mkdirs()) {
            throw new IOException("Unable to write to directory " + value);
        }
        return value;
    }

    @NotNull
    public static String nonOptionalFile(@NotNull CommandLine cmd, @NotNull String param)   {
        return fileExists(nonOptionalValue(cmd, param));
    }

    @NotNull
    public static String fileExists(@NotNull String fullPath) {
        boolean exists = Files.exists(new File(fullPath).toPath());
        if (!exists) {
            throw new IllegalArgumentException(format("Path [%s] did not exist and is mandatory.", fullPath));
        }
        return fullPath;
    }
}
