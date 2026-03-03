package com.hartwig.hmftools.finding.datamodel;

import org.jspecify.annotations.Nullable;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record Virus(
        @NotNull DriverFields driver,
        @NotNull String name,
        @NotNull VirusBreakendQCStatus qcStatus,
        int integrations,
        @Nullable OncogenicVirus oncogenicVirus,
        double percentageCovered,
        double meanCoverage,
        @Nullable Double expectedClonalCoverage
) implements Driver
{
    public enum VirusBreakendQCStatus
    {
        NO_ABNORMALITIES,
        LOW_VIRAL_COVERAGE,
        EXCESSIVE_VIRAL_COVERAGE,
        ASSEMBLY_DOWNSAMPLED,
        CHILD_TAXID_REFERENCE,
        UNCLEAR_TAXID_ASSIGNMENT
    }

    public enum OncogenicVirus
    {
        MCV,
        EBV,
        HPV,
        HBV,
        HHV8
    }

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

    @Override
    public double driverLikelihood()
    {
        return driver.driverLikelihood();
    }
}
