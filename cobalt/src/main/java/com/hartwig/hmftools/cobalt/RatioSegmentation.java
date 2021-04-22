package com.hartwig.hmftools.cobalt;

import java.io.IOException;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.utils.pcf.PCFFile;
import com.hartwig.hmftools.common.utils.r.RExecutor;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class RatioSegmentation {

    private static final Logger LOGGER = LogManager.getLogger(RatioSegmentation.class);

    private final String outputDirectory;
    private final ExecutorService executorService;

    RatioSegmentation(final ExecutorService executorService, final String outputDirectory) {
        this.outputDirectory = outputDirectory;
        this.executorService = executorService;
    }

    void applySegmentation(@NotNull final String reference, @NotNull final String tumor) throws ExecutionException, InterruptedException {
        final String ratioFile = CobaltRatioFile.generateFilenameForReading(outputDirectory, tumor);
        final List<Future<Object>> futures = Lists.newArrayList();
        futures.add(executorService.submit(() -> ratioSegmentation(ratioFile, reference, "referenceGCDiploidRatio")));
        futures.add(executorService.submit(() -> ratioSegmentation(ratioFile, tumor, "tumorGCRatio")));

        for (Future<Object> future : futures) {
            future.get();
        }

        LOGGER.info("Segmentation Complete");
    }

    @Nullable
    private Object ratioSegmentation(@NotNull final String ratioFile, @NotNull final String sample, @NotNull final String column)
            throws IOException, InterruptedException {
        final String pcfFile = PCFFile.generateRatioFilename(outputDirectory, sample);
        int result = RExecutor.executeFromClasspath("r/ratioSegmentation.R", ratioFile, column, pcfFile);
        if (result != 0) {
            throw new IOException("R execution failed. Unable to complete segmentation.");
        }

        return null;
    }
}
