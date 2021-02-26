package com.hartwig.hmftools.common.lims;

import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class LimsChecker {

    private static final Logger LOGGER = LogManager.getLogger(LimsChecker.class);

    public static boolean checkViralInsertions(@Nullable LimsJsonSampleData sampleData, @Nullable LimsCohortConfig cohort,
            @NotNull String sampleId) {
        if (sampleData != null && cohort != null) {
            if (sampleData.reportViralInsertions()) {
                if (!cohort.reportViral()) {
                    LOGGER.warn("Consent of viral insertions is true, but must be false for sample '{}'", sampleId);
                }
                return true;
            } else {
                if (cohort.reportViral()) {
                    LOGGER.warn("Consent of viral insertions is false, but must be true for sample '{}'", sampleId);
                }
                return false;
            }
        } else {
            return false;
        }
    }

    public static boolean checkGermlineVariants(@Nullable LimsJsonSampleData sampleData, @Nullable LimsCohortConfig cohort,
            @NotNull String sampleId) {
        if (sampleData != null && cohort != null) {
            if (sampleData.reportGermlineVariants()) {
                if (!cohort.reportGermline()) {
                    LOGGER.warn("Consent of report germline variants is true, but must be false for sample '{}'", sampleId);
                }
                return true;
            } else {
                if (cohort.reportGermline()) {
                    LOGGER.warn("Consent of report germline variants is false, but must be true for sample '{}'", sampleId);
                }
                return false;
            }
        } else {
            return false;
        }
    }
}
