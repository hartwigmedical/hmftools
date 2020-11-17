package com.hartwig.hmftools.serve.sources.iclusion;

import java.util.List;

import com.hartwig.hmftools.iclusion.data.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.serve.sources.ExtractionOutput;

import org.jetbrains.annotations.NotNull;

public class IclusionExtractor {

    @NotNull
    public ExtractionOutput extractFromIclusionTrials(@NotNull List<IclusionTrial> trials) {

        for (IclusionTrial trial : trials) {
            for (IclusionMutationCondition mutationCondition : trial.mutationConditions()) {
                mutationCondition.logicType();
            }
        }
        return null;
    }
}
