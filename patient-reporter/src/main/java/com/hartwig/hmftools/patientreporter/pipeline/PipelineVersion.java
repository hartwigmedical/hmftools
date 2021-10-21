package com.hartwig.hmftools.patientreporter.pipeline;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.patientreporter.algo.AnalysedPatientReporter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PipelineVersion {

    private static final Logger LOGGER = LogManager.getLogger(PipelineVersion.class);

    @VisibleForTesting
    public static void checkPipelineVersion(@Nullable String actualPipelineVersion, @NotNull String expectedPipelineVersion,
            boolean overridePipelineVersion) {
        if (overridePipelineVersion) {
            if (actualPipelineVersion == null) {
                LOGGER.warn("No known pipeline version is known!");
            }
            LOGGER.warn("Pipeline version is overridden! The version is {} and the expected version is {}",
                    actualPipelineVersion,
                    expectedPipelineVersion);
        } else {
            if (actualPipelineVersion != null && !actualPipelineVersion.equals(expectedPipelineVersion)) {
                throw new IllegalArgumentException(
                        "The expected pipeline version " + expectedPipelineVersion + " is different from the actual pipeline version "
                                + actualPipelineVersion + "!");
            } else if (actualPipelineVersion == null) {
                throw new IllegalArgumentException(
                        "No pipeline version is known! The expected pipeline version is " + expectedPipelineVersion);
            }
        }
    }
}
