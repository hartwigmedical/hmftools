package com.hartwig.hmftools.amber;

import java.io.IOException;

import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.pcf.PCFFile;
import com.hartwig.hmftools.common.r.RExecutor;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class BAFSegmentation {

    private static final Logger LOGGER = LogManager.getLogger(BAFSegmentation.class);

    @NotNull
    private final String outputDirectory;

    BAFSegmentation(@NotNull final String outputDirectory) {
        this.outputDirectory = outputDirectory;
    }

    void applySegmentation(@NotNull final String tumor) throws InterruptedException, IOException {
        final String ratioFile = AmberBAFFile.generateAmberFilename(outputDirectory, tumor);
        bafSegmentation(ratioFile, tumor);
        LOGGER.info("BAFSegmentation Complete");
    }

    private void bafSegmentation(@NotNull final String ratioFile, @NotNull final String sample) throws IOException, InterruptedException {
        final String pcfFile = PCFFile.generateBAFFilename(outputDirectory, sample);
        int result = RExecutor.executeFromClasspath("r/bafSegmentation.R", ratioFile, pcfFile);
        if (result != 0) {
            throw new IOException("R execution failed. Unable to complete segmentation.");
        }
    }
}
