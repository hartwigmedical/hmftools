package com.hartwig.hmftools.protect.evidence;

import com.hartwig.hmftools.common.sv.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.variant.DriverInterpretation;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class EvidenceDriverLikelihood {

    @Nullable
    public static Boolean interpretVariant(@NotNull DriverInterpretation driverInterpretation) {
        switch (driverInterpretation) {
            case HIGH:
                return true;
            case MEDIUM:
            case LOW:
                return false;
            default:
                return null;
        }
    }

    @Nullable
    public static Boolean interpretSignatures() {
        return null;
    }

    @Nullable
    public static Boolean interpretFusions(@NotNull FusionLikelihoodType likelihood) {
        switch (likelihood) {
            case HIGH:
                return true;
            case LOW:
                return false;
            case NA:
                return null;
            default:
                return null;
        }
    }

    public static boolean interpretDisruptions() {
        // All homozygous disruption has a high driver likelihood
        return true;
    }

    public static boolean interpretCNV() {
        //All loss/gains has a high driver likelihood
        return true;
    }

    @Nullable
    public static Boolean interpretChord() {
        return null;
    }

    public static boolean interpretVirus() {
        //We use only high driver viruses for evidence matching
        return true;
    }
}