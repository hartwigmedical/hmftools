package com.hartwig.hmftools.orange.algo.protect;

import java.util.List;

import com.hartwig.hmftools.common.protect.ProtectEvidence;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class EvidenceEvaluator {

    private EvidenceEvaluator() {
    }

    public static boolean hasEvidence(@NotNull List<ProtectEvidence> evidences, @Nullable String gene, @NotNull String event) {
        for (ProtectEvidence evidence : evidences) {
            if ((gene == null || gene.equals(evidence.gene())) && evidence.event().equals(event)) {
                return true;
            }
        }

        return false;
    }
}
