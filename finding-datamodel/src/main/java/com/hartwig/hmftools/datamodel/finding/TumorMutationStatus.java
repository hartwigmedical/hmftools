package com.hartwig.hmftools.datamodel.finding;

import com.hartwig.hmftools.datamodel.purple.PurpleTumorMutationalStatus;

import org.jetbrains.annotations.NotNull;

import io.soabase.recordbuilder.core.RecordBuilder;

@RecordBuilder
public record TumorMutationStatus(
        @NotNull String findingKey,
        double tumorMutationalBurdenPerMb,
        @NotNull PurpleTumorMutationalStatus tumorMutationalBurdenStatus,
        int tumorMutationalLoad,
        @NotNull PurpleTumorMutationalStatus tumorMutationalLoadStatus,
        int svTumorMutationalBurden
) implements Finding {}
