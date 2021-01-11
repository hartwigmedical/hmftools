package com.hartwig.hmftools.protect.evidence;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidenceItem;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.jetbrains.annotations.NotNull;

public final class ProtectEvidenceItemTestFactory {

    private ProtectEvidenceItemTestFactory() {
    }

    @NotNull
    public static ImmutableProtectEvidenceItem.Builder createDefault(boolean onLabel, @NotNull EvidenceDirection direction,
            @NotNull EvidenceLevel level) {
        return ImmutableProtectEvidenceItem.builder()
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
