package com.hartwig.hmftools.cobalt;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.hartwig.hmftools.common.utils.pcf.PCFFile;
import com.hartwig.hmftools.common.utils.r.RExecutor;

public class RatioSegmentation
{
    public static void applyRatioSegmentation(
            final ExecutorService executorService, final String outputDir, final String ratioFile,
            final String reference, final String tumor, int gamma) throws ExecutionException, InterruptedException
    {
        final List<Future<Object>> futures = new ArrayList<>();

        if(reference != null)
        {
            futures.add(executorService.submit(() -> ratioSegmentation(outputDir, ratioFile, reference, "referenceGCDiploidRatio", gamma)));
        }
        if(tumor != null)
        {
            futures.add(executorService.submit(() -> ratioSegmentation(outputDir, ratioFile, tumor, "tumorGCRatio", gamma)));
        }

        for(Future<Object> future : futures)
        {
            future.get();
        }

        CB_LOGGER.info("Segmentation Complete");
    }

    private static Object ratioSegmentation(
            final String outputDir, final String ratioFile, final String sample, final String column, int gamma)
            throws IOException, InterruptedException
    {
        final String pcfFile = PCFFile.generateRatioFilename(outputDir, sample);
        int result = RExecutor.executeFromClasspath("r/ratioSegmentation.R", ratioFile, column, pcfFile, String.valueOf(gamma));
        if(result != 0)
        {
            throw new IOException("R execution failed. Unable to complete segmentation.");
        }

        return null;
    }
}
