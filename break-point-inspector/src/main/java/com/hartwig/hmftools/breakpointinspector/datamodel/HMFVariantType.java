package com.hartwig.hmftools.breakpointinspector.datamodel;

import org.jetbrains.annotations.NotNull;

public enum HMFVariantType {
    INS,
    DEL,
    INV3,
    INV5,
    DUP;

    @NotNull
    public static String orientation(@NotNull HMFVariantType type) {
        switch (type) {
            case INS:
            case DEL:
                return "INNIE";
            case INV3:
                return "TANDEM_RIGHT";
            case INV5:
                return "TANDEM_LEFT";
            case DUP:
                return "OUTIE";
            default:
                return "ERROR";
        }
    }
}
