package com.hartwig.hmftools.protect.evidence;

import com.hartwig.hmftools.common.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.variant.DriverInterpretation;

import org.jetbrains.annotations.NotNull;

public final class EvidenceDriverLikelihood {

    private EvidenceDriverLikelihood() {
    }

    public static boolean interpretVariant(@NotNull DriverInterpretation driverInterpretation) {
        switch (driverInterpretation) {
            case HIGH:
                return true;
            case MEDIUM:
            case LOW:
                return false;
            default:
                throw new IllegalStateException("Unrecognized driver variant interpretation: " + driverInterpretation);
        }
    }

    public static boolean interpretFusion(@NotNull FusionLikelihoodType likelihood) {
        switch (likelihood) {
            case HIGH:
                return true;
            case LOW:
            case NA:
                return false;
            default:
                throw new IllegalStateException("Unrecognized fusion likelihood type: " + likelihood);
        }
    }

    public static boolean interpretCopyNumber() {
        // All loss/gains have a high driver likelihood
        return true;
    }

    public static boolean interpretVirus() {
        // We use only high driver viruses for evidence matching
        return true;
    }

    public static boolean interpretWildType() {
       // All wild types genes has a low driver likelihood
        return false;
    }
}