package com.hartwig.hmftools.serve.sources.iclusion;

import java.util.List;

import com.hartwig.hmftools.iclusion.data.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.data.IclusionMutationLogicType;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.serve.sources.ExtractionOutput;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class IclusionExtractor {

    private static final Logger LOGGER = LogManager.getLogger(IclusionExtractor.class);

    @NotNull
    public ExtractionOutput extractFromIclusionTrials(@NotNull List<IclusionTrial> trials) {
        // We assume filtered trials.
        for (IclusionTrial trial : trials) {
            assert !trial.acronym().isEmpty() && !trial.mutationConditions().isEmpty();
            for (IclusionMutationCondition mutationCondition : trial.mutationConditions()) {
                assert mutationCondition.logicType() == IclusionMutationLogicType.OR;
                // TODO implement.
            }
        }
        return null;
    }
}
