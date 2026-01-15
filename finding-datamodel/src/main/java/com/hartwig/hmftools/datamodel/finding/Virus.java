package com.hartwig.hmftools.datamodel.finding;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.datamodel.virus.VirusInterpretation;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import io.soabase.recordbuilder.core.RecordBuilder;

@RecordBuilder
public record Virus(
        @NotNull DriverFields driver,
        @NotNull String name,
        @NotNull VirusBreakendQCStatus qcStatus,
        int integrations,
        @Nullable VirusInterpretation interpretation,
        double percentageCovered,
        double meanCoverage,
        @Nullable Double expectedClonalCoverage
) implements Driver
{
    @NotNull
    @Override
    public String findingKey()
    {
        return driver.findingKey();
    }

    @NotNull
    @Override
    public DriverSource driverSource()
    {
        return driver.driverSource();
    }

    @NotNull
    @Override
    public ReportedStatus reportedStatus()
    {
        return driver.reportedStatus();
    }

    @NotNull
    @Override
    public DriverInterpretation driverInterpretation()
    {
        return driver.driverInterpretation();
    }
}
