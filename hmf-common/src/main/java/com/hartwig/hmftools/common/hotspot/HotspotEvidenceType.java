package com.hartwig.hmftools.common.hotspot;

import org.jetbrains.annotations.NotNull;

public enum HotspotEvidenceType {
    INFRAME,
    KNOWN,
    SNV,
    MNV,
    INSERT,
    DELETE;

    @NotNull
    static HotspotEvidenceType fromVariantHotspot(@NotNull VariantHotspot hotspot) {
        if (hotspot.isSimpleDelete()) {
            return DELETE;
        }

        if (hotspot.isSimpleInsert()) {
            return INSERT;
        }

        if (hotspot.isMNV()) {
            return MNV;
        }

        return SNV;
    }

}
