package com.hartwig.hmftools.common.hotspot;

import org.jetbrains.annotations.NotNull;

public enum HotspotEvidenceType {
    INFRAME,
    SNV,
    MNV,
    INSERT,
    DELETE;

    @NotNull
    static HotspotEvidenceType fromVariantHotspot(@NotNull VariantHotspot hotspot) {
        if (hotspot.isDelete()) {
            return DELETE;
        }

        if (hotspot.isInsert()) {
            return INSERT;
        }

        if (hotspot.isMNV()) {
            return MNV;
        }

        return SNV;
    }

}
