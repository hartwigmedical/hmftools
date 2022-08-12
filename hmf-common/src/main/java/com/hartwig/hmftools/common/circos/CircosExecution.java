package com.hartwig.hmftools.common.circos;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.util.CollectionUtil;

public class CircosExecution
{
    private static final Logger LOGGER = LogManager.getLogger(CircosExecution.class);

    private final String executable;

    public CircosExecution(@NotNull final String executable)
    {
        this.executable = executable;
    }

    @Nullable
    public Integer generateCircos(@NotNull final String inputConfig, @NotNull final String outputPath, @NotNull final String outputFile,
            @NotNull final String errorPath) throws IOException, InterruptedException
    {
        final File redirectErrorFile = new File(errorPath + File.separator + outputFile + ".error");
        final File redirectOutputFile = new File(errorPath + File.separator + outputFile + ".out");

        final String[] command = new String[8];
        command[0] = executable;
        command[1] = "-nosvg";
        command[2] = "-conf";
        command[3] = new File(inputConfig).getAbsolutePath();
        command[4] = "-outputdir";
        command[5] = new File(outputPath).getAbsolutePath();
        command[6] = "-outputfile";
        command[7] = outputFile;

        LOGGER.info(String.format("Generating " + outputFile + " via command: %s", CollectionUtil.join(Arrays.asList(command), " ")));
        int result = new ProcessBuilder(command).redirectError(redirectErrorFile).redirectOutput(redirectOutputFile).start().waitFor();
        if(result != 0)
        {
            LOGGER.error("Fatal error creating circos plot. Examine error file " + redirectErrorFile.toString() + " for details.");
            return 0;
        }

        final File finalFile = new File(outputPath + File.separator + outputFile);
        if(!finalFile.exists())
        {
            LOGGER.error("Failed to create file {}", finalFile.toString());
        }

        return 0;
    }
}
