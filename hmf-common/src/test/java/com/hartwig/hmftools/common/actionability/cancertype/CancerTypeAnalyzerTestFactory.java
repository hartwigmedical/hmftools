package com.hartwig.hmftools.common.actionability.cancertype;

import java.io.IOException;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class CancerTypeAnalyzerTestFactory {

    private CancerTypeAnalyzerTestFactory() {
    }

    @NotNull
    public static CancerTypeAnalyzer buildWithOneCancerTypeMapping(@NotNull String cancerType, @NotNull String doid) {
        CancerTypeToDOIDMappingEntry cancerTypeToDOIDMappingEntry =
                ImmutableCancerTypeToDOIDMappingEntry.builder().cancerType(cancerType).addDoids(doid).build();

        PrimaryTumorToDOIDMapping primaryTumorToDOIDMapping;
        try {
            primaryTumorToDOIDMapping = PrimaryTumorToDOIDMapping.createFromResource();
        } catch (IOException exception) {
            throw new IllegalStateException("Should always be able to create production mapping");
        }

        return new CancerTypeAnalyzer(Lists.newArrayList(cancerTypeToDOIDMappingEntry), primaryTumorToDOIDMapping);
    }
}
