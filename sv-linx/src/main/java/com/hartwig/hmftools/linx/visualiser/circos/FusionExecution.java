package com.hartwig.hmftools.linx.visualiser.circos;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.utils.r.RExecutor;
import com.hartwig.hmftools.linx.visualiser.CircosConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class FusionExecution
{
    private static final Logger LOGGER = LogManager.getLogger(FusionExecution.class);

    private final String plotFile;
    private final String exonFile;
    private final String proteinDomainFile;

    public FusionExecution(final String sampleName, final String imageName, final String dataDir, final String plotDir)
    {
        this.plotFile = plotDir + File.separator + imageName;
        this.exonFile = dataDir + File.separator + sampleName + ".fusions.tsv";
        this.proteinDomainFile = dataDir + File.separator + sampleName + ".protein_domains.tsv";
    }

    public Integer executeR(CircosConfig config, double labelSize) throws IOException, InterruptedException
    {
        int result = RExecutor.executeFromClasspath("r/fusionPlot.R",
                proteinDomainFile,
                exonFile,
                plotFile,
                String.valueOf(labelSize),
                String.valueOf(config.FusionLegendRows),
                String.valueOf(config.FusionLegendHeightPerRow),
                String.valueOf(config.FusionHeight)
        );
        if (result != 0)
        {
            LOGGER.warn("Error adding fusion plots.");
        }

        return result;
    }

}
