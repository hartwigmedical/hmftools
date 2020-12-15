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
    public static ImmutableProtectEvidenceItem.Builder createDefault(boolean onLabel, EvidenceDirection direction, EvidenceLevel level) {
        return ImmutableProtectEvidenceItem.builder()
                .source(Knowledgebase.VICC_CGI)
                .direction(direction)
                .level(level)
                .treatment("treatment")
                .onLabel(onLabel)
                .genomicEvent("event")
                .reported(true)
                .urls(Sets.newHashSet("url"));
    }
}
