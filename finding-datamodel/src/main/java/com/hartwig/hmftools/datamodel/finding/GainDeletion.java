package com.hartwig.hmftools.datamodel.finding;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record GainDeletion(
        @NotNull DriverFields driver,
        @NotNull Type type,
        @NotNull CopyNumberInterpretation interpretation,
        @NotNull String gene,
        @NotNull String chromosome,
        @NotNull String chromosomeBand,
        @NotNull String transcript,
        boolean isCanonical,
        double tumorMinCopies,
        double tumorMaxCopies,
        double tumorMinMinorAlleleCopies,
        double chromosomeArmCopies) implements Driver
{
    public enum Type
    {
        GERMLINE_DEL_HOM_IN_TUMOR,
        GERMLINE_DEL_HET_IN_TUMOR,
        GERMLINE_GAIN,
        SOMATIC_GAIN,
        SOMATIC_DEL,
        SOMATIC_LOH
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
