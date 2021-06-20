package com.hartwig.hmftools.linx.visualiser.circos;

import static com.hartwig.hmftools.linx.visualiser.SvVisualiser.VIS_LOGGER;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import com.hartwig.hmftools.common.utils.r.RExecutor;
import com.hartwig.hmftools.linx.visualiser.CircosConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.util.IOUtils;
import org.jetbrains.annotations.NotNull;

public class ChromosomeRangeExecution
{
    private final String mPlotFile;
    private final String mBandFile;
    private final String mChromosomeFile;

    public ChromosomeRangeExecution(final String sampleName, final String imageName, final String dataDir, final String plotDir)
    {
        mPlotFile = plotDir + File.separator + imageName;
        mBandFile = dataDir + File.separator + sampleName + ".cytoBand.txt";
        mChromosomeFile = dataDir + File.separator + sampleName + ".chromosome.circos";
    }

    public Integer executeR(@NotNull final CircosConfig config, double labelSize) throws IOException, InterruptedException
    {
        writeCytobands();
        int result = RExecutor.executeFromClasspath("r/chromosomeRangePlot.R",
                mChromosomeFile, mBandFile, mPlotFile, String.valueOf(labelSize),
                String.valueOf(config.ChromosomeRangeHeight), String.valueOf(config.ChromosomeRangeColumns));

        if (result != 0)
        {
            VIS_LOGGER.warn("Error adding chromosomal context");
        }

        return result;
    }

    private void writeCytobands() throws IOException
    {
        final String template = readResource("/r/cytoBand.txt");
        Files.write(new File(mBandFile).toPath(), template.getBytes(StandardCharsets.UTF_8));
    }

    @NotNull
    private String readResource(@NotNull final String resource) throws IOException
    {
        InputStream in = getClass().getResourceAsStream(resource);
        BufferedReader reader = new BufferedReader(new InputStreamReader(in));
        return IOUtils.toString(reader);
    }

}
