package com.hartwig.hmftools.protect.evidence;

import com.hartwig.hmftools.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.sources.Source;

import org.jetbrains.annotations.NotNull;

public class ProtectEvidenceItemTest {

    @NotNull
    public static ImmutableProtectEvidenceItem.Builder createDefault(boolean onLabel, EvidenceDirection direction, EvidenceLevel level) {
        return ImmutableProtectEvidenceItem.builder()
                .source(Source.CGI)
                .direction(direction)
                .level(level)
                .treatment("treatment")
                .onLabel(onLabel)
                .genomicEvent("event")
                .reported(true)
                .url("url");
    }

}
