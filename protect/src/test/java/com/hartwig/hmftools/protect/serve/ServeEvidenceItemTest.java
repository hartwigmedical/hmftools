package com.hartwig.hmftools.protect.serve;

import com.hartwig.hmftools.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.sources.Source;

import org.jetbrains.annotations.NotNull;

public class ServeEvidenceItemTest {


    @NotNull
    public static ImmutableServeEvidenceItem.Builder createDefault(boolean onLabel, EvidenceDirection direction, EvidenceLevel level) {
        return ImmutableServeEvidenceItem.builder()
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
