package com.hartwig.hmftools.protect.evidence;

import com.hartwig.hmftools.common.protect.ImmutableProtectEvidenceItem;
import com.hartwig.hmftools.common.serve.EvidenceDirection;
import com.hartwig.hmftools.common.serve.EvidenceLevel;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.jetbrains.annotations.NotNull;

public class ProtectEvidenceItemTest {

    @NotNull
    public static ImmutableProtectEvidenceItem.Builder createDefault(boolean onLabel, EvidenceDirection direction, EvidenceLevel level) {
        return ImmutableProtectEvidenceItem.builder()
                .source(Knowledgebase.CGI)
                .direction(direction)
                .level(level)
                .treatment("treatment")
                .onLabel(onLabel)
                .genomicEvent("event")
                .reported(true)
                .url("url");
    }

}
