package com.hartwig.hmftools.linx.visualiser.circos;

import static com.hartwig.hmftools.linx.visualiser.SvVisualiser.VIS_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import com.hartwig.hmftools.common.genome.chromosome.CytoBands;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.RExecutor;
import com.hartwig.hmftools.linx.visualiser.CircosConfig;

public class ChromosomeRangeExecution
{
    private final String mPlotFile;
    private final String mBandFile;
    private final String mChromosomeFile;
    private final RefGenomeVersion mRefGenomeVersion;

    public ChromosomeRangeExecution(
            final String sampleName, final RefGenomeVersion rgVersion, final String imageName, final String dataDir, final String plotDir)
    {
        mPlotFile = plotDir + File.separator + imageName;
        mBandFile = dataDir + File.separator + sampleName + ".cytoBand.txt";
        mChromosomeFile = dataDir + File.separator + sampleName + ".chromosome.circos";
        mRefGenomeVersion = rgVersion;
    }

    public Integer executeR(final CircosConfig config, double labelSize) throws IOException, InterruptedException
    {
        writeCytobands();

        int result = RExecutor.executeFromClasspath("r/chromosomeRangePlot.R",
                mChromosomeFile, mBandFile, mPlotFile, String.valueOf(labelSize),
                String.valueOf(config.ChromosomeRangeHeight), String.valueOf(config.ChromosomeRangeColumns));

        if(result != 0)
        {
            VIS_LOGGER.warn("error adding chromosomal context");
        }

        return result;
    }

    private void writeCytobands() throws IOException
    {
        final String template = CytoBands.resourceAsString(mRefGenomeVersion);
        Files.write(new File(mBandFile).toPath(), template.getBytes(StandardCharsets.UTF_8));
    }
}
