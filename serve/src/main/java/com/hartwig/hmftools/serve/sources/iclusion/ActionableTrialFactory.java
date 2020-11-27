package com.hartwig.hmftools.serve.sources.iclusion;

import java.util.List;

import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.iclusion.data.IclusionTumorLocation;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ActionableTrialFactory {

    private ActionableTrialFactory() {
    }

    @NotNull
    public static List<ActionableTrial> toActionableTrials(@NotNull IclusionTrial trial) {
        // TODO URL should be something like "https://iclusion.org/hmf/" + ext"
        ImmutableActionableTrial.Builder actionableBuilder = ImmutableActionableTrial.builder()
                .source(Knowledgebase.ICLUSION)
                .treatment(trial.acronym())
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .url(Strings.EMPTY);

        // TODO Make sure DOIDs are present. Currently not part of iClusion API
        List<ActionableTrial> actionableTrials = Lists.newArrayList();
        for (IclusionTumorLocation tumorLocation : trial.tumorLocations()) {
            actionableTrials.add(actionableBuilder.cancerType(tumorLocation.primaryTumorLocation()).doid(Strings.EMPTY).build());
        }

        return actionableTrials;
    }
}
