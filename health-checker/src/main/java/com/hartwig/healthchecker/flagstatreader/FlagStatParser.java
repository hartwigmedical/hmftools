package com.hartwig.healthchecker.flagstatreader;

import java.io.IOException;

import com.hartwig.healthchecker.exception.EmptyFileException;

import org.jetbrains.annotations.NotNull;

interface FlagStatParser {

    @NotNull
    FlagStatData parse(@NotNull String filePath, @NotNull String filter)
                    throws IOException, EmptyFileException;
}
