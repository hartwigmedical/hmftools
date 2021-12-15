package com.hartwig.hmftools.serve.sources.actin;

import java.util.List;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class ActionableTrialFactory {

    @NotNull
    public List<ActionableTrial> toActionableEntries(@NotNull ActinTrial actionTrial) {
        List<ActionableTrial> trialEntries = Lists.newArrayList();

        com.hartwig.hmftools.serve.sources.actin.ImmutableActionableTrial.Builder actionableBuilder =
                com.hartwig.hmftools.serve.sources.actin.ImmutableActionableTrial.builder()
                        .source(Knowledgebase.ACTIN)
                        .treatment(actionTrial.trialId() + "(" + actionTrial.cohortId() + ")")
                        .level(EvidenceLevel.B)
                        .direction(EvidenceDirection.RESPONSIVE)
                        .urls(Sets.newHashSet(""));

        trialEntries.add(actionableBuilder.cancerType("").doid("").build());

        return trialEntries;
    }
}
