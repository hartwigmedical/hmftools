package com.hartwig.hmftools.common.protect;

import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ProtectTestFactory {

    private ProtectTestFactory() {
    }

    @NotNull
    public static ProtectEvidence createTestProtectEvidence() {
        return testEvidenceBuilder().build();
    }

    @NotNull
    public static ImmutableProtectEvidence.Builder testEvidenceBuilder() {
        return ImmutableProtectEvidence.builder()
                .genomicEvent(Strings.EMPTY)
                .germline(false)
                .reported(true)
                .treatment(Strings.EMPTY)
                .onLabel(false)
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE);
    }
}
