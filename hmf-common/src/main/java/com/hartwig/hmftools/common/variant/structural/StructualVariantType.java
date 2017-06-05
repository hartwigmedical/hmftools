package com.hartwig.hmftools.common.variant.structural;

public enum StructualVariantType {

    BND, INV, DEL, INS, DUP;

    public static StructualVariantType fromAttribute(String svType) {
        if (svType.startsWith("DUP")) {
            return DUP;
        }
        return StructualVariantType.valueOf(svType);
    }
}
