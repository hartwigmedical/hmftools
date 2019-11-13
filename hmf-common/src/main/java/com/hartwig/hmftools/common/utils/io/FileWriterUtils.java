package com.hartwig.hmftools.common.utils.io;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

import org.jetbrains.annotations.NotNull;

public final class FileWriterUtils {

    private FileWriterUtils() {
    }

    @NotNull
    public static BufferedWriter createBufferedWriter(final String outputFile, boolean appendIfExists) throws IOException {
        Path outputPath = Paths.get(outputFile);

        if (Files.exists(outputPath))
        {
            if (appendIfExists)
            {
                return Files.newBufferedWriter(outputPath, StandardOpenOption.APPEND);
            }
            else
            {
                return Files.newBufferedWriter(outputPath, StandardOpenOption.TRUNCATE_EXISTING);
            }
        }
        else
        {
            return Files.newBufferedWriter(outputPath, StandardOpenOption.CREATE);
        }
    }

    public static void closeBufferedWriter(BufferedWriter writer)
    {
        if(writer == null)
            return;

        try
        {
            writer.close();
        }
        catch (IOException e)
        {
            throw new IllegalStateException("Could not close buffered writer: " + writer + ": " + e.getMessage());
        }
    }
}
