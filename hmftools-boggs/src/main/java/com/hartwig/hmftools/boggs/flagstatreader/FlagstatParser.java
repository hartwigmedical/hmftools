package com.hartwig.hmftools.boggs.flagstatreader;

import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;

public interface FlagstatParser {

    @NotNull
    FlagstatData parse(@NotNull File file) throws IOException;
}
