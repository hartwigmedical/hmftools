package com.hartwig.hmftools.iclusion.qc;

import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class IclusionTrialChecker {

    private static final Logger LOGGER = LogManager.getLogger(IclusionTrialChecker.class);

    private IclusionTrialChecker() {
    }

    public static void check(@NotNull Iterable<IclusionTrial> trials) {
        LOGGER.info("Performing QC on iClusion trials");

        for (IclusionTrial trial : trials) {
            if (trial.acronym().isEmpty()) {
                LOGGER.warn("Empty acronym for trial with id {} and title '{}'", trial.id(), trial.title());
            }

            if (!trial.eudra().isEmpty() && (!trial.eudra().contains("-") || trial.eudra().contains(" "))) {
                // EUDRA codes are formatted as 'yyyy-xxxxxx-xx' where 'yyyy' is (full) year, eg 2019.
                LOGGER.warn("Potentially incorrect EUDRA code found for {} with id {}: '{}'", trial.acronym(), trial.id(), trial.eudra());
            }
        }
    }
}
