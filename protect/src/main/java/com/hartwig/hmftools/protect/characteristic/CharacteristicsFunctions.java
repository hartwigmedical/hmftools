package com.hartwig.hmftools.protect.characteristic;

import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicsComparator;

import org.jetbrains.annotations.NotNull;

public final class CharacteristicsFunctions {

    private CharacteristicsFunctions() {
    }

    public static boolean hasExplicitCutoff(@NotNull ActionableCharacteristic signature) {
        if (signature.comparator() == null) {
            return false;
        }
        else if (signature.comparator() == TumorCharacteristicsComparator.BETWEEN) {
            return signature.minCutoff() != null && signature.maxCutoff() != null;
        } else {
            return signature.maxCutoff() != null;
        }
    }

    public static boolean evaluateVersusCutoff(@NotNull ActionableCharacteristic signature, double value) {
        assert hasExplicitCutoff(signature);

        switch (signature.comparator()) {
            case EQUAL_OR_LOWER:
                return value <= signature.maxCutoff();
            case LOWER:
                return value < signature.maxCutoff();
            case EQUAL_OR_GREATER:
                return value >= signature.maxCutoff();
            case GREATER:
                return value > signature.maxCutoff();
            case BETWEEN:
                return value >= signature.minCutoff() && value <= signature.maxCutoff();
            default: {
                throw new IllegalStateException("Unrecognized comparator: " + signature.comparator());
            }
        }
    }
}