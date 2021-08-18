package com.hartwig.hmftools.orange.algo.selection;

import java.util.List;

import com.hartwig.hmftools.common.protect.ProtectEvidence;

import org.jetbrains.annotations.NotNull;

public final class EvidenceSelector {

    private EvidenceSelector() {
    }

    public static boolean hasEvidence(@NotNull List<ProtectEvidence> evidences, @NotNull String genomicEvent) {
        for (ProtectEvidence evidence : evidences) {
            if (evidence.genomicEvent().equals(genomicEvent)) {
                return true;
            }
        }

        return false;
    }
}
