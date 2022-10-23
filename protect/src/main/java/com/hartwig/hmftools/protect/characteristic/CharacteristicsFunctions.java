package com.hartwig.hmftools.protect.characteristic;

import com.hartwig.serve.datamodel.characteristic.ActionableCharacteristic;

import org.jetbrains.annotations.NotNull;

public final class CharacteristicsFunctions {

    private CharacteristicsFunctions() {
    }

    public static boolean hasExplicitCutoff(@NotNull ActionableCharacteristic signature) {
        return signature.comparator() != null && signature.cutoff() != null;
    }

    public static boolean evaluateVersusCutoff(@NotNull ActionableCharacteristic signature, double value) {
        assert hasExplicitCutoff(signature);

        switch (signature.comparator()) {
            case EQUAL_OR_LOWER:
                return value <= signature.cutoff();
            case LOWER:
                return value < signature.cutoff();
            case EQUAL_OR_GREATER:
                return value >= signature.cutoff();
            case GREATER:
                return value > signature.cutoff();
            default: {
                throw new IllegalStateException("Unrecognized comparator: " + signature.comparator());
            }
        }
    }
}