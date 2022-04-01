package com.hartwig.hmftools.common.protect;

import com.beust.jcommander.internal.Sets;
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
                .event(Strings.EMPTY)
                .evidenceType(ProtectEvidenceType.ANY_MUTATION)
                .germline(false)
                .reported(true)
                .treatment(Strings.EMPTY)
                .onLabel(false)
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Sets.newHashSet())
                .protectSources(ImmutableProtectSource.builder().build());
    }
}
