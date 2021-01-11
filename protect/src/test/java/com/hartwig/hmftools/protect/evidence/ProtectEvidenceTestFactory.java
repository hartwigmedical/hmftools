package com.hartwig.hmftools.protect.evidence;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.actionability.ActionabilityTestUtil;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ProtectEvidenceTestFactory {

    private ProtectEvidenceTestFactory() {
    }

    @NotNull
    public static ActionableEvent createTestBaseEvent() {
        return ActionabilityTestUtil.create(Knowledgebase.VICC_CGI,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                EvidenceLevel.A,
                EvidenceDirection.RESPONSIVE,
                Sets.newHashSet());
    }

    @NotNull
    public static ImmutableProtectEvidence.Builder createDefault(boolean onLabel, @NotNull EvidenceDirection direction,
            @NotNull EvidenceLevel level) {
        return ImmutableProtectEvidence.builder()
                .genomicEvent("event")
                .germline(false)
                .source(Knowledgebase.VICC_CGI)
                .reported(true)
                .treatment("treatment")
                .onLabel(onLabel)
                .level(level)
                .direction(direction)
                .urls(Sets.newHashSet("url"));
    }
}
