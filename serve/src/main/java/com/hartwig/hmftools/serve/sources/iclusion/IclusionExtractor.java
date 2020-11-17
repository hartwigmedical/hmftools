package com.hartwig.hmftools.serve.sources.iclusion;

import java.util.List;

import com.hartwig.hmftools.iclusion.data.IclusionMutation;
import com.hartwig.hmftools.iclusion.data.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.serve.sources.ExtractionOutput;
import com.hartwig.hmftools.serve.sources.ImmutableExtractionOutput;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class IclusionExtractor {

    private static final Logger LOGGER = LogManager.getLogger(IclusionExtractor.class);

    @NotNull
    public ExtractionOutput extractFromIclusionTrials(@NotNull List<IclusionTrial> trials) {
        // We assume filtered trials (no empty acronyms, only OR mutations, and no negated mutations

        for (IclusionTrial trial : trials) {
            for (IclusionMutationCondition mutationCondition : trial.mutationConditions()) {
                for (IclusionMutation mutation : mutationCondition.mutations()) {
                    LOGGER.info("Interpreting '{}' on '{}' for {}", mutation.name(), mutation.gene(), trial.acronym());
                }
            }
        }
        return ImmutableExtractionOutput.builder().build();
    }
}
