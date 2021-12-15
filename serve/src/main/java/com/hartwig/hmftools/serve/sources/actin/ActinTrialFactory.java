package com.hartwig.hmftools.serve.sources.actin;

import java.util.List;

import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class ActinTrialFactory {

    @NotNull
    public List<ActinTrial> toActionableEntries(@NotNull ActinEntry actionTrial) {
        List<ActinTrial> trialEntries = Lists.newArrayList();

        ImmutableActinTrial.Builder actionableBuilder = ImmutableActinTrial.builder()
                .source(Knowledgebase.ACTIN)
                .treatment(actionTrial.trial())
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE);

        trialEntries.add(actionableBuilder.cancerType("").doid("").build());

        return trialEntries;
    }
}
