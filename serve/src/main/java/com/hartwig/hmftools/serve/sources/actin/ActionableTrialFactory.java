package com.hartwig.hmftools.serve.sources.actin;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.jetbrains.annotations.NotNull;

public class ActionableTrialFactory {

    @NotNull
    public Set<ActionableTrial> toActionableEntries(@NotNull ActinTrial actionTrial) {
        Set<ActionableTrial> trialEntries = Sets.newHashSet();

        com.hartwig.hmftools.serve.sources.actin.ImmutableActionableTrial.Builder actionableBuilder =
                com.hartwig.hmftools.serve.sources.actin.ImmutableActionableTrial.builder()
                        .source(Knowledgebase.ICLUSION)
                        .treatment("")
                        .level(EvidenceLevel.B)
                        .direction(EvidenceDirection.RESPONSIVE)
                        .urls(Sets.newHashSet(""));

        trialEntries.add(actionableBuilder.cancerType("").doid("").build());

        return trialEntries;
    }
}
