package com.hartwig.hmftools.puritypatho.variants;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.io.reader.FileReader;

import org.jetbrains.annotations.NotNull;

public class ReadingFileVariantDetection {

    @NotNull
    public static List<String> read(@NotNull final String fileName) throws IOException {
        return FileReader.build().readLines(new File(fileName).toPath());
    }
}
