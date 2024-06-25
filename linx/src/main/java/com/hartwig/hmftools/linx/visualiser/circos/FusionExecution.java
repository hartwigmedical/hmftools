package com.hartwig.hmftools.linx.visualiser.circos;

import static com.hartwig.hmftools.linx.visualiser.SvVisualiser.VIS_LOGGER;
import static com.hartwig.hmftools.linx.visualiser.circos.FusionDataWriter.FUSION_PLOT_TSV;
import static com.hartwig.hmftools.linx.visualiser.circos.FusionDataWriter.PROTEIN_PLOT_TSV;

import java.io.IOException;

import com.hartwig.hmftools.common.utils.r.RExecutor;
import com.hartwig.hmftools.linx.visualiser.CircosConfig;

public class FusionExecution
{
    private final String mPlotFile;
    private final String mExonFile;
    private final String mProteinDomainFile;

    public FusionExecution(final String sampleName, final String imageName, final String dataDir, final String plotDir)
    {
        mPlotFile = plotDir + imageName;
        mExonFile = dataDir + sampleName + FUSION_PLOT_TSV;
        mProteinDomainFile = dataDir + sampleName + PROTEIN_PLOT_TSV;
    }

    public Integer executeR(CircosConfig config, double labelSize) throws IOException, InterruptedException
    {
        int result = RExecutor.executeFromClasspath("r/fusionPlot.R",
                mProteinDomainFile, mExonFile, mPlotFile, String.valueOf(labelSize), String.valueOf(config.FusionLegendRows),
                String.valueOf(config.FusionLegendHeightPerRow), String.valueOf(config.FusionHeight));

        if(result != 0)
        {
            VIS_LOGGER.error("error adding fusion plots");
        }

        return result;
    }

}
