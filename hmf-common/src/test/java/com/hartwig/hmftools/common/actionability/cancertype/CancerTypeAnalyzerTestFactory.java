package com.hartwig.hmftools.common.actionability.cancertype;

import java.io.IOException;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public final class CancerTypeAnalyzerTestFactory {

    private CancerTypeAnalyzerTestFactory() {
    }

    @NotNull
    public static CancerTypeAnalyzer buildWithOneCancerTypeMapping(@NotNull String cancerType, @NotNull String doid) {
        PrimaryTumorToDOIDMapper primaryTumorToDOIDMapper;
        try {
            primaryTumorToDOIDMapper = PrimaryTumorToDOIDMapper.createFromResource();
        } catch (IOException exception) {
            throw new IllegalStateException("Should always be able to create production mapping");
        }

        Map<String, Set<String>> doidsPerCancerType = Maps.newHashMap();
        doidsPerCancerType.put(cancerType, Sets.newHashSet(doid));

        return new CancerTypeAnalyzer(new CancerTypeToDOIDMapper(doidsPerCancerType), primaryTumorToDOIDMapper);
    }
}
