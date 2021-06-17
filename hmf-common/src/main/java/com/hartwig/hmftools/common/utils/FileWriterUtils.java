package com.hartwig.hmftools.common.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Map;

import com.google.common.collect.Maps;

import org.apache.commons.cli.CommandLine;
import org.jetbrains.annotations.NotNull;

public final class FileWriterUtils
{
    public static final String OUTPUT_DIR = "output_dir";

    public static String parseOutputDir(@NotNull final CommandLine cmd)
    {
        String outputDir = cmd.getOptionValue(OUTPUT_DIR);
        if(outputDir == null)
            return null;

        return checkAddDirSeparator(outputDir);
    }

    public static boolean checkCreateOutputDir(final String outputDirPath)
    {
        if (Files.exists(Paths.get(outputDirPath)))
            return true;

        final File outputDir = new File(outputDirPath);
        return outputDir.mkdirs();
    }

    public static String checkAddDirSeparator(@NotNull final String outputDir)
    {
        if(outputDir.endsWith(File.separator))
            return outputDir;

        return outputDir + File.separator;
    }

    @NotNull
    public static BufferedWriter createBufferedWriter(final String outputFile, boolean appendIfExists) throws IOException
    {
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

    public static Map<String,Integer> createFieldsIndexMap(final String fieldsHeader, final String delimiter)
    {
        final String[] items = fieldsHeader.split(delimiter,-1);
        final Map<String,Integer> fieldsIndexMap = Maps.newLinkedHashMap();

        for(int i = 0; i < items.length; ++i)
        {
            fieldsIndexMap.put(items[i], i);
        }

        return fieldsIndexMap;
    }

}
