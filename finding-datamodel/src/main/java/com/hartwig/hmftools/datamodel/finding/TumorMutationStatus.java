package com.hartwig.hmftools.datamodel.finding;

import com.hartwig.hmftools.datamodel.purple.PurpleTumorMutationalStatus;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record TumorMutationStatus(
        @NotNull String findingKey,
        double tumorMutationalBurdenPerMb,
        @NotNull PurpleTumorMutationalStatus tumorMutationalBurdenStatus,
        int tumorMutationalLoad,
        @NotNull PurpleTumorMutationalStatus tumorMutationalLoadStatus,
        int svTumorMutationalBurden
) implements Finding
{
}
