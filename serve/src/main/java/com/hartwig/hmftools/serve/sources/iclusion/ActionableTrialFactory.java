package com.hartwig.hmftools.serve.sources.iclusion;

import java.util.List;

import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTumorLocation;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public final class ActionableTrialFactory {

    private ActionableTrialFactory() {
    }

    @NotNull
    public static List<ActionableTrial> toActionableTrials(@NotNull IclusionTrial trial) {
        ImmutableActionableTrial.Builder actionableBuilder = ImmutableActionableTrial.builder()
                .source(Knowledgebase.ICLUSION)
                .treatment(trial.acronym())
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .url("https://iclusion.org/hmf/" + trial.id());

        List<ActionableTrial> actionableTrials = Lists.newArrayList();
        for (IclusionTumorLocation tumorLocation : trial.tumorLocations()) {
            for (String doid : tumorLocation.doids()) {
                actionableTrials.add(actionableBuilder.cancerType(tumorLocation.primaryTumorLocation()).doid(doid).build());
            }
        }

        return actionableTrials;
    }
}
