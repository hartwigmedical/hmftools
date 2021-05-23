package com.hartwig.hmftools.amber;

import java.io.IOException;

import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.utils.pcf.PCFFile;
import com.hartwig.hmftools.common.utils.r.RExecutor;

import org.jetbrains.annotations.NotNull;

public class BAFSegmentation
{
    private final String mOutputDir;

    public BAFSegmentation(final String outputDir)
    {
        mOutputDir = outputDir;
    }

    public void applySegmentation(final String tumor) throws InterruptedException, IOException
    {
        final String ratioFile = AmberBAFFile.generateAmberFilenameForReading(mOutputDir, tumor);
        final String pcfFile = PCFFile.generateBAFFilename(mOutputDir, tumor);
        int result = RExecutor.executeFromClasspath("r/bafSegmentation.R", ratioFile, pcfFile);
        if(result != 0)
        {
            throw new IOException("R execution failed. Unable to complete segmentation.");
        }
    }
}
