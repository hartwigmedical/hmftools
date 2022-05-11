package com.hartwig.hmftools.protect.evidence;

import com.hartwig.hmftools.common.sv.linx.FusionLikelihoodType;
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

    public static boolean interpretFusions(@NotNull FusionLikelihoodType likelihood) {
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

    public static boolean interpretCNV() {
        //All loss/gains has a high driver likelihood
        return true;
    }

    public static boolean interpretVirus() {
        //We use only high driver viruses for evidence matching
        return true;
    }
}