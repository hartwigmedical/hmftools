package com.hartwig.hmftools.patientdb.clinical.lims;

import com.hartwig.hmftools.patientdb.clinical.lims.cohort.LimsCohortConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class LimsChecker {

    private static final Logger LOGGER = LogManager.getLogger(LimsChecker.class);

    private LimsChecker() {
    }

    public static boolean checkViralInsertions(@Nullable LimsJsonSampleData sampleData, @Nullable LimsCohortConfig cohort,
            @NotNull String sampleId) {
        if (sampleData != null && cohort != null) {
            if (sampleData.reportViralPresence()) {
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

    @Nullable
    public static String toHospitalPathologySampleIdForReport(@NotNull String hospitalPathologySampleId, @NotNull String tumorSampleId,
            @NotNull LimsCohortConfig cohortConfig) {
        if (cohortConfig.requireHospitalPAId()) {
            if (!hospitalPathologySampleId.equals(Lims.NOT_AVAILABLE_STRING) && !hospitalPathologySampleId.isEmpty()
                    && isValidHospitalPathologySampleId(hospitalPathologySampleId)) {
                return hospitalPathologySampleId;
            } else {

                LOGGER.warn("Missing or invalid hospital pathology sample ID for sample '{}': {}. Please fix!",
                        tumorSampleId,
                        hospitalPathologySampleId);

                return null;
            }
        } else {
            if (!hospitalPathologySampleId.isEmpty() && !hospitalPathologySampleId.equals(Lims.NOT_AVAILABLE_STRING)) {
                LOGGER.info("Skipping hospital pathology sample ID for sample '{}': {}", hospitalPathologySampleId, tumorSampleId);
            }

            return null;
        }
    }

    public static void checkHospitalPatientId(@NotNull String hospitalPatientId, @NotNull String sampleId,
            @NotNull LimsCohortConfig cohortConfig) {
        if (cohortConfig.requireHospitalId()) {
            if (hospitalPatientId.equals(Lims.NOT_AVAILABLE_STRING) || hospitalPatientId.isEmpty()) {
                LOGGER.warn("Missing hospital patient sample ID for sample '{}': {}. Please fix!", sampleId, hospitalPatientId);
            }
        }
    }

    private static boolean isValidHospitalPathologySampleId(@NotNull String hospitalPathologySampleId) {
        boolean tMatch;
        boolean cMatch;

        if (hospitalPathologySampleId.split("-")[1].length() <= 6 && hospitalPathologySampleId.startsWith("T") && !hospitalPathologySampleId
                .startsWith("C")) {
            tMatch = hospitalPathologySampleId.substring(1, 3).matches("[0-9]+") && hospitalPathologySampleId.charAt(3) == '-'
                    && hospitalPathologySampleId.substring(4, 4 + hospitalPathologySampleId.split("-")[1].length()).matches("[0-9]+");
        } else if (hospitalPathologySampleId.contains(" ")) {
            tMatch = hospitalPathologySampleId.substring(1, 3).matches("[0-9]+") && hospitalPathologySampleId.charAt(3) == '-'
                    && hospitalPathologySampleId.substring(4, 4 + hospitalPathologySampleId.split("-")[1].split("\\s")[0].length())
                    .matches("[0-9]+") && hospitalPathologySampleId.split("\\s")[1].contains("I");
        } else {
            tMatch = false;
        }

        if (hospitalPathologySampleId.split("-")[1].length() <= 6 && hospitalPathologySampleId.startsWith("C") && !hospitalPathologySampleId
                .startsWith("T")) {
            cMatch = hospitalPathologySampleId.substring(1, 3).matches("[0-9]+") && hospitalPathologySampleId.charAt(3) == '-'
                    && hospitalPathologySampleId.substring(4, 4 + hospitalPathologySampleId.split("-")[1].length()).matches("[0-9]+");
        } else if (hospitalPathologySampleId.contains(" ")) {
            cMatch = hospitalPathologySampleId.substring(1, 3).matches("[0-9]+") && hospitalPathologySampleId.charAt(3) == '-'
                    && hospitalPathologySampleId.substring(4, 4 + hospitalPathologySampleId.split("-")[1].split("\\s")[0].length())
                    .matches("[0-9]+");
        } else {
            cMatch = false;
        }

        return tMatch || cMatch;
    }
}
