package com.hartwig.hmftools.orange.algo.protect;

import java.util.List;

import com.hartwig.hmftools.common.protect.ProtectEvidence;

import org.jetbrains.annotations.NotNull;

public final class ProtectInterpreter {

    private ProtectInterpreter() {
    }

    @NotNull
    public static ProtectInterpretedData interpret(@NotNull List<ProtectEvidence> evidences) {
        return ImmutableProtectInterpretedData.builder()
                .reportableEvidences(EvidenceSelector.reportableEvidence(evidences))
                .reportableTrials(EvidenceSelector.reportableTrials(evidences))
                .unreportedEvidences(EvidenceSelector.unreportedEvidence(evidences))
                .unreportedTrials(EvidenceSelector.unreportedTrials(evidences))
                .build();
    }
}
