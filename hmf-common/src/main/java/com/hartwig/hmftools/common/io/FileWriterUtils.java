package com.hartwig.hmftools.common.io;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

public class FileWriterUtils
{
    public static BufferedWriter createBufferedWriter(final String outputFile, boolean appendIfExists)
    {
        try
        {
            Path outputPath = Paths.get(outputFile);

            if(Files.exists(outputPath))
            {
                if(appendIfExists)
                    return Files.newBufferedWriter(outputPath, StandardOpenOption.APPEND);
                else
                    return Files.newBufferedWriter(outputPath, StandardOpenOption.TRUNCATE_EXISTING);
            }
            else
            {
                return Files.newBufferedWriter(outputPath, StandardOpenOption.CREATE);
            }

        } catch (IOException e)
        {
            return null;
        }
    }

    public static void closeBufferedWriter(final BufferedWriter writer)
    {
        if(writer == null)
            return;

        try
        {
            writer.close();
        }
        catch (IOException e)
        {

        }
    }
}
