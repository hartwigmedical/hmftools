package com.hartwig.hmftools.finding.datamodel;

import com.hartwig.hmftools.finding.datamodel.driver.Driver;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFields;
import com.hartwig.hmftools.finding.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.finding.datamodel.driver.DriverSource;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record ChromosomeArmCopyNumber(
        @NotNull DriverFields driver,
        @NotNull String chromosome,
        @NotNull ChromosomeArm arm,
        @NotNull Type type,
        double copyNumber,
        double relativeCopyNumber
) implements Driver
{
    public enum ChromosomeArm
    {
        P,
        Q,
    }

    public enum Type
    {
        GAIN,
        LOSS,
        DIPLOID
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

    @Override
    public boolean isReported()
    {
        return driver.isReported();
    }
}
