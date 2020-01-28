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
                // CCMO code should always be present (dutch trial code).
                LOGGER.warn("No CCMO code configured for trial with acronym '{}'", trial.acronym());
            }

            if (!trial.eudra().isEmpty() && (!trial.eudra().contains("-") || trial.eudra().contains(" "))) {
                // EUDRA codes are formatted as 'yyyy-xxxxxx-xx' where 'yyyy' is (full) year, eg 2019.
                LOGGER.warn("Potentially incorrect EUDRA code found for {}: '{}'", trial.acronym(), trial.eudra());
            }
         }
    }
}
