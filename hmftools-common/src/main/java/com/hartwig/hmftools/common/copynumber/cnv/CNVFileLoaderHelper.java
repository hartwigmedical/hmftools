package com.hartwig.hmftools.common.copynumber.cnv;

import java.io.File;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CNVFileLoaderHelper {

    // KODU: copynumber data is stored in {run}/copyNumber/{sampleR}_{sampleT}/freec/{sampleT}<>.bam_CNVs
    private static final String COPYNUMBER_BASE_DIRECTORY = "copyNumber";
    private static final String COPYNUMBER_SAMPLE_CONNECTOR = "_";
    private static final String FREEC_ALGO_DIRECTORY = "freec";

    private CNVFileLoaderHelper() {
    }

    @NotNull
    public static String getFreecBasePath(@NotNull final String runDirectory, @NotNull final String refSample,
            @Nullable final String tumorSample) {
        final String baseDir = runDirectory + File.separator + COPYNUMBER_BASE_DIRECTORY + File.separator;
        final String sampleDir =
                tumorSample != null ? refSample + COPYNUMBER_SAMPLE_CONNECTOR + tumorSample : refSample;

        return baseDir + sampleDir + File.separator + FREEC_ALGO_DIRECTORY;
    }
}
