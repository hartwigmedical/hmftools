package com.hartwig.hmftools.patientreporter.algo;

import java.util.Comparator;

import com.hartwig.hmftools.common.protect.ProtectEvidence;

import org.jetbrains.annotations.NotNull;

public class ProtectComparator implements Comparator<ProtectEvidence>  {

    @Override
    public int compare(@NotNull ProtectEvidence protectEvidence1, @NotNull ProtectEvidence protectEvidence2) {
        if (protectEvidence1.treatment().equals(protectEvidence2.treatment())) {
            if (protectEvidence1.level().equals(protectEvidence2.level())) {
                return protectEvidence1.compareTo(protectEvidence2);
            }

        }

        return protectEvidence2.compareTo(protectEvidence1);
    }
}
