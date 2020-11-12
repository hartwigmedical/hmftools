package com.hartwig.hmftools.common.clinical;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.lims.LimsStudy;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class PatientPrimaryTumorFunctions {

    private static final Logger LOGGER = LogManager.getLogger(PatientPrimaryTumorFunctions.class);

    private PatientPrimaryTumorFunctions() {
    }

    @Nullable
    public static PatientPrimaryTumor findPrimaryTumorForSample(@NotNull List<PatientPrimaryTumor> patientPrimaryTumors,
            @NotNull String sample) {
        String patientIdentifier = toPatientIdentifier(sample);

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

    @NotNull
    private static String toPatientIdentifier(@NotNull String sample) {
        LimsStudy study = LimsStudy.fromSampleId(sample);
        if (sample.length() >= 12 && (study == LimsStudy.CPCT || study == LimsStudy.DRUP || study == LimsStudy.WIDE
                || study == LimsStudy.CORE)) {
            return sample.substring(0, 12);
        } else if (sample.toUpperCase().startsWith("COLO829")) {
            return "COLO829";
        }

        // If we want to fetch the primary tumor for a sample we cannot derive patient of, we assume sample and patient are equal.
        return sample;
    }
}
