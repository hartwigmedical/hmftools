package com.hartwig.hmftools.cobalt.segment;

import java.io.IOException;

import com.hartwig.hmftools.common.purple.pcf.PCFFile;
import com.hartwig.hmftools.common.purple.ratio.ReadRatioFile;
import com.hartwig.hmftools.common.r.RExecutor;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class PCFSegment {

    private static final Logger LOGGER = LogManager.getLogger(PCFSegment.class);

    private final String outputDirectory;

    public PCFSegment(final String outputDirectory) {
        this.outputDirectory = outputDirectory;
    }

    public void ratioSegmentation(@NotNull final String sample) throws IOException, InterruptedException {
        final String ratioFile = ReadRatioFile.generateFilename(outputDirectory, sample);
        final String pcfFile = PCFFile.generateRatioFilename(outputDirectory, sample);
        int result = RExecutor.executeFromClasspath("r/ratioSegmentation.R", ratioFile, pcfFile);
        if (result != 0) {
            throw new IOException("R execution failed. Unable to complete segmentation.");
        }

        LOGGER.info("Segmentation Complete");
    }

}
