package com.hartwig.hmftools.protect.actionability.cancertype;

import java.io.IOException;
import java.util.Collections;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CancerTypeAnalyzer {

    @NotNull
    private final CancerTypeToDOIDMapper cancerTypeToDOIDMapper;
    @NotNull
    private final PrimaryTumorToDOIDMapper primaryTumorToDOIDMapper;

    @NotNull
    public static CancerTypeAnalyzer createFromKnowledgeBase(@NotNull String knowledgebaseCancerTypesTsv) throws IOException {
        return new CancerTypeAnalyzer(CancerTypeToDOIDMapper.createFromFile(knowledgebaseCancerTypesTsv),
                PrimaryTumorToDOIDMapper.createFromResource());
    }

    @VisibleForTesting
    CancerTypeAnalyzer(@NotNull final CancerTypeToDOIDMapper cancerTypeToDOIDMapper,
            @NotNull final PrimaryTumorToDOIDMapper primaryTumorToDOIDMapper) {
        this.cancerTypeToDOIDMapper = cancerTypeToDOIDMapper;
        this.primaryTumorToDOIDMapper = primaryTumorToDOIDMapper;
    }

    public boolean isCancerTypeMatch(@NotNull String knowledgebaseCancerType, @Nullable String primaryTumorLocation) {
        if (primaryTumorLocation == null) {
            return false;
        }

        Set<String> doidsForPrimaryTumorLocation = primaryTumorToDOIDMapper.findDoids(primaryTumorLocation);
        if (doidsForPrimaryTumorLocation == null) {
            return false;
        }

        Set<String> doidsForCancerType = cancerTypeToDOIDMapper.findDoids(knowledgebaseCancerType);
        return doidsForCancerType != null && !Collections.disjoint(doidsForPrimaryTumorLocation, doidsForCancerType);
    }
}
