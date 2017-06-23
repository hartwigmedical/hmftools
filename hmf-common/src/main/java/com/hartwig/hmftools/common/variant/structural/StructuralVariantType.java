package com.hartwig.hmftools.common.variant.structural;

public enum StructuralVariantType {

    BND, INV, DEL, INS, DUP;

    public static StructuralVariantType fromAttribute(String svType) {
        if (svType.startsWith("DUP")) {
            return DUP;
        }
        return StructuralVariantType.valueOf(svType);
    }
}
