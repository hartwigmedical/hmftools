package com.hartwig.hmftools.boggs.flagstatreader;

import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;

public interface FlagstatParser2 {

    @NotNull
    FlagstatData2 parse(@NotNull File file) throws IOException;
}
