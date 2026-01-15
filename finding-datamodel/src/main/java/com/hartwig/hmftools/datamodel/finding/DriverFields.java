package com.hartwig.hmftools.datamodel.finding;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;

import org.jetbrains.annotations.NotNull;

import io.soabase.recordbuilder.core.RecordBuilder;

@RecordBuilder
public record DriverFields(
        @NotNull String findingKey,
        @NotNull DriverSource driverSource,
        @NotNull ReportedStatus reportedStatus,
        @NotNull DriverInterpretation driverInterpretation
) {}
