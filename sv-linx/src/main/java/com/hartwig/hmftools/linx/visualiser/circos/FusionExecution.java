package com.hartwig.hmftools.linx.visualiser.circos;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.r.RExecutor;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class FusionExecution
{
    private static final Logger LOGGER = LogManager.getLogger(FusionExecution.class);

    private final String plotFile;
    private final String exonFile;
    private final String proteinDomainFile;

    public FusionExecution(final String sample, final String dataDir, final String plotDir)
    {
        this.plotFile = plotDir + File.separator + sample + ".png";
        this.exonFile = dataDir + File.separator + sample + ".fusions.tsv";
        this.proteinDomainFile = dataDir + File.separator + sample + ".protein_domains.tsv";
    }

    public Integer executeR() throws IOException, InterruptedException
    {
        int result = RExecutor.executeFromClasspath("r/fusionPlot.R",
                proteinDomainFile,
                exonFile,
                plotFile);
        if (result != 0)
        {
            LOGGER.warn("Error generating R copy number plots.");
        }

        return result;
    }

}
