package com.hartwig.hmftools.common.circos;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

public class CircosExecution
{
    private static final Logger LOGGER = LogManager.getLogger(CircosExecution.class);

    private final String executable;

    public CircosExecution(final String executable)
    {
        this.executable = executable;
    }

    @Nullable
    public Integer generateCircos(
            final String inputConfig, final String outputPath, final String outputFile,
            final String errorPath) throws IOException, InterruptedException
    {
        final File redirectOutputFile = new File(errorPath + File.separator + outputFile + ".out");

        final List<String> command = List.of(
            executable,
            "-nosvg",
            "-conf",
            new File(inputConfig).getAbsolutePath(),
            "-outputdir",
            new File(outputPath).getAbsolutePath(),
            "-outputfile",
            outputFile);

        LOGGER.info(String.format("Generating " + outputFile + " via command: %s", String.join(" ", command)));
        Process process = new ProcessBuilder(command).redirectOutput(redirectOutputFile).start();
        int result = process.waitFor();

        if(result != 0)
        {
            System.err.print(new String(process.getErrorStream().readAllBytes()));
            LOGGER.error("Fatal error creating circos plot.");
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
