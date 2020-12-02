package com.hartwig.hmftools.iclusion.classification;

import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.EventClassifierFactory;
import com.hartwig.hmftools.common.serve.classification.MutationType;
import com.hartwig.hmftools.iclusion.datamodel.IclusionMutation;

import org.jetbrains.annotations.NotNull;

public final class MutationTypeExtractor {

    @NotNull
    private static final EventClassifier CLASSIFIER = EventClassifierFactory.buildClassifier(IclusionClassificationConfig.build());

    private MutationTypeExtractor() {
    }

    @NotNull
    public static MutationType classify(@NotNull IclusionMutation mutation) {
        return CLASSIFIER.determineType(mutation.gene(), mutation.name());
    }
}
