package com.hartwig.hmftools.iclusion.classification;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.serve.classification.EventMatcher;
import com.hartwig.hmftools.common.serve.classification.MutationType;
import com.hartwig.hmftools.iclusion.data.IclusionMutation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class MutationTypeClassifier {

    private static final Logger LOGGER = LogManager.getLogger(MutationTypeClassifier.class);

    @NotNull
    private static final Map<MutationType, EventMatcher> MATCHERS = Maps.newHashMap();

    public MutationTypeClassifier() {
    }

    @NotNull
    public static MutationType classify(@NotNull IclusionMutation mutation) {
        return MutationType.UNKNOWN;
    }
}
