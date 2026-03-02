package com.hartwig.hmftools.finding.datamodel;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;

import org.jspecify.annotations.Nullable;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record GainDeletion(
        @NotNull DriverFields driver,
        @NotNull String gene,
        @NotNull String chromosome,
        @NotNull String chromosomeBand,
        @NotNull String transcript,
        boolean isCanonical,
        @Nullable Type germlineType,
        @NotNull Type somaticType,
        @NotNull GeneExtent geneExtent,
        @Nullable ExonRange exonRange, // null if exon range info not available
        double tumorMinCopies,
        double tumorMaxCopies,
        double tumorMinMinorAlleleCopies,
        double chromosomeArmCopies,
        @Nullable Double germlineMinCopyNumber) implements Driver
{
    public enum Type
    {
        GAIN,
        DEL,
        HEL_DEL,
        CN_NEUTRAL_LOH
    }

    public enum GeneExtent
    {
        WHOLE_GENE,
        PARTIAL_GENE
    }

    public record ExonRange(
            @Nullable Integer exonStart, // null if starts before first exon
            @Nullable Integer exonEnd    // null if ends after last exon
    )
    {}

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
