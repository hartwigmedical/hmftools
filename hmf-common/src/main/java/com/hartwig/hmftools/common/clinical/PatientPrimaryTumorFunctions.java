package com.hartwig.hmftools.common.clinical;

import java.util.List;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class PatientPrimaryTumorFunctions {

    private static final Logger LOGGER = LogManager.getLogger(PatientPrimaryTumorFunctions.class);

    private PatientPrimaryTumorFunctions() {
    }

    @Nullable
    public static PatientPrimaryTumor findPrimaryTumorForPatient(@NotNull List<PatientPrimaryTumor> patientPrimaryTumors,
            @NotNull String patientIdentifier) {
        LOGGER.info(patientIdentifier);
        List<PatientPrimaryTumor> matchingIdPrimaryTumors = patientPrimaryTumors.stream()
                .filter(patientPrimaryTumor -> patientPrimaryTumor.patientIdentifier().equals(patientIdentifier))
                .collect(Collectors.toList());

        // We should never have more than one primary tumor for a single patient.
        assert matchingIdPrimaryTumors.size() < 2;

        if (matchingIdPrimaryTumors.size() == 1) {
            return matchingIdPrimaryTumors.get(0);
        } else {
            LOGGER.warn("Could not find patient '{}' in list of primary tumors!", patientIdentifier);
            return null;
        }
    }
}
