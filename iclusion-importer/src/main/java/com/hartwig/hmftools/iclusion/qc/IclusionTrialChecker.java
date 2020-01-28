package com.hartwig.hmftools.iclusion.qc;

import com.hartwig.hmftools.iclusion.data.IclusionTrial;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class IclusionTrialChecker {

    private static final Logger LOGGER = LogManager.getLogger(IclusionTrialChecker.class);

    private IclusionTrialChecker() {
    }

    public static void check(@NotNull Iterable<IclusionTrial> trials) {
        for (IclusionTrial trial : trials) {
            if (trial.acronym().isEmpty()) {
                LOGGER.warn("Empty acronym for trial with title '{}'", trial.title());
            }

            if (trial.ccmo().isEmpty()) {
                LOGGER.warn("No CCMO code configured for trial with acronym '{}'", trial.acronym());
            }
        }
    }
}
