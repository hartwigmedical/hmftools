package com.hartwig.hmftools.common.circos;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;
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
            final String inputConfig, final String outputPath, final String outputFile) throws IOException, InterruptedException
    {
        String plotFilePath = outputPath + File.separator + outputFile;

        // we have to delete existing plot file first, otherwise circos could silently fail
        File plotFile = new File(plotFilePath);
        if(plotFile.exists())
        {
            plotFile.delete();
        }

        final List<String> command = List.of(
            executable,
            "-nosvg",
            "-conf",
            new File(inputConfig).getAbsolutePath(),
            "-outputdir",
            new File(outputPath).getAbsolutePath(),
            "-outputfile",
            outputFile);

        LOGGER.info(String.format("generating " + outputFile + " via command: %s", String.join(" ", command)));

        ProcessBuilder processBuilder = new ProcessBuilder(command);
        Map<String, String> environment = processBuilder.environment();

        // Circos writes the temporary file circos.colorlist. However, this file would be overwritten multiple times when invoking circos
        // from multiple threads which results in a race condition. We make sure that circos writes to different temporary directories to
        // avoid this race condition
        Path tempDir = Files.createTempDirectory("circos");
        environment.put("TMPDIR", tempDir.toString());

        // must redirect error stream to stdout, as circos print some errors to stdout
        Process process = processBuilder.redirectErrorStream(true).start();

        int result = process.waitFor();

        FileUtils.deleteDirectory(tempDir.toFile());

        if(result != 0)
        {
            System.err.print(new String(process.getInputStream().readAllBytes()));
            LOGGER.error("Fatal error creating circos plot.");
            return 0;
        }

        plotFile = new File(outputPath + File.separator + outputFile);
        if(!plotFile.exists())
        {
            System.err.print(new String(process.getInputStream().readAllBytes()));
            LOGGER.error("Failed to create file {}", plotFile);
        }

        return 0;
    }
}
