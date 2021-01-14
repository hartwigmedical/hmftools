package com.hartwig.hmftools.patientdb.interpretator.drup;

import org.jetbrains.annotations.Nullable;

public class BaselineInterpretator {

    private BaselineInterpretator() {

    }

    @Nullable
    public static String interpretationBaselinePrimaryTumorLocation(@Nullable String primaryTumorLocationBastType,
            @Nullable String primaryTumorLocationBastTypeOther, @Nullable String primaryTumorReg) {
        String primaryTumorLocation = null;
        if (primaryTumorLocationBastType != null && !primaryTumorLocationBastType.isEmpty()) {
            if (primaryTumorLocationBastType.equals("Other, specify")) {
                primaryTumorLocation = primaryTumorLocationBastTypeOther;
            } else {
                primaryTumorLocation = primaryTumorLocationBastType;
            }
        }

        // See DEV-1713 for why below choices have been made.
        if (primaryTumorReg != null && !primaryTumorReg.isEmpty()) {
            if (primaryTumorLocation != null) {
                primaryTumorLocation = primaryTumorLocation + " + " + primaryTumorReg;
            } else {
                primaryTumorLocation = primaryTumorReg;
            }
        }
        return primaryTumorLocation;
    }

    @Nullable
    public static String setFinalPrimaryTumorLocation(@Nullable String primaryTumorCohort, @Nullable String primaryTumorLocation) {
        // See DEV-1713 for why below choices have been made.
        String finalPrimaryTumorLocation;
        if (primaryTumorCohort != null && !primaryTumorCohort.isEmpty()) {
            String lowerPrimaryTumorCohort = primaryTumorCohort.trim().toLowerCase();
            if (primaryTumorLocation != null && (lowerPrimaryTumorCohort.contains("biliary tract") || lowerPrimaryTumorCohort.contains(
                    "colon") || lowerPrimaryTumorCohort.contains("urinary organ") || lowerPrimaryTumorCohort.contains("head, face and neck")
                    || lowerPrimaryTumorCohort.contains("salivary gland"))) {
                finalPrimaryTumorLocation = primaryTumorCohort + " + " + primaryTumorLocation;
            } else {
                finalPrimaryTumorLocation = primaryTumorCohort;
            }
        } else {
            finalPrimaryTumorLocation = primaryTumorLocation;
        }
        return finalPrimaryTumorLocation;
    }
}
