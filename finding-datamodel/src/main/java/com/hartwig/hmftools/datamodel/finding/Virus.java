package com.hartwig.hmftools.datamodel.finding;

import java.util.Set;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;
import com.hartwig.hmftools.datamodel.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.datamodel.virus.VirusInterpretation;

import org.jspecify.annotations.Nullable;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

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
    public String event()
    {
        return driver.event();
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

    @Override
    public double driverLikelihood()
    {
        return driver.driverLikelihood();
    }

    @Override
    public boolean isReported()
    {
        return driver.isReported();
    }

    @Override
    public Set<String> genes()
    {
        return Set.of();
    }
}
