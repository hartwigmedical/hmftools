package com.hartwig.hmftools.linx.visualiser.circos;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import com.hartwig.hmftools.common.utils.r.RExecutor;
import com.hartwig.hmftools.linx.visualiser.SvCircosConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.util.IOUtils;
import org.jetbrains.annotations.NotNull;

public class ChromosomeRangeExecution
{
    private static final Logger LOGGER = LogManager.getLogger(ChromosomeRangeExecution.class);

    private final String plotFile;
    private final String bandFile;
    private final String chromosomeFile;

    public ChromosomeRangeExecution(final String sampleName, final String imageName, final String dataDir, final String plotDir)
    {
        this.plotFile = plotDir + File.separator + imageName;
        this.bandFile = dataDir + File.separator + sampleName + ".cytoBand.txt";
        this.chromosomeFile = dataDir + File.separator + sampleName + ".chromosome.circos";
    }

    public Integer executeR(@NotNull final SvCircosConfig config, double labelSize) throws IOException, InterruptedException
    {
        writeCytobands();
        int result = RExecutor.executeFromClasspath("r/chromosomeRangePlot.R",
                chromosomeFile,
                bandFile,
                plotFile,
                String.valueOf(labelSize),
                String.valueOf(config.chromosomeRangeHeight()),
                String.valueOf(config.chromosomeRangeColumns())
        );
        if (result != 0)
        {
            LOGGER.warn("Error adding chromosomal context");
        }

        return result;
    }

    private void writeCytobands() throws IOException
    {
        final String template = readResource("/r/cytoBand.txt");
        Files.write(new File(bandFile).toPath(), template.getBytes(StandardCharsets.UTF_8));
    }

    @NotNull
    private String readResource(@NotNull final String resource) throws IOException
    {
        InputStream in = getClass().getResourceAsStream(resource);
        BufferedReader reader = new BufferedReader(new InputStreamReader(in));
        return IOUtils.toString(reader);
    }

}
