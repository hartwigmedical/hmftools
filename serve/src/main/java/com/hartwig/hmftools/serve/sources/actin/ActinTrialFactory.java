package com.hartwig.hmftools.serve.sources.actin;

import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;

import org.jetbrains.annotations.NotNull;

public final class ActinTrialFactory {

    private ActinTrialFactory() {
    }

    @NotNull
    public static ActinTrial toActinTrial(@NotNull ActinEntry actionTrial) {
        return ImmutableActinTrial.builder()
                .source(Knowledgebase.ACTIN)
                .treatment(actionTrial.trial())
                .level(EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE)
                .cancerType("Advanced Solid Tumor")
                .doid("162")
                .build();
    }
}
