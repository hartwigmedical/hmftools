package com.hartwig.hmftools.finding.datamodel;

import org.jspecify.annotations.Nullable;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

@RecordBuilder
@RecordBuilder.Options(nullableAnnotationClass = "org.jspecify.annotations.Nullable", defaultNotNull = true)
public record GainDeletion(
        @NotNull DriverFields driver,
        @NotNull String gene,
        @NotNull String chromosome,
        @NotNull String chromosomeBand,
        @NotNull String transcript,
        boolean isCanonical,
        @NotNull Type somaticType,
        @Nullable Type germlineType,
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
        HOM_DEL,
        HET_DEL,
        CN_NEUTRAL_LOH,
        NONE
    }

    public enum GeneExtent
    {
        FULL_GENE,
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

    public boolean isLossOfHeterozygosity()
    {
        return (somaticType() == Type.CN_NEUTRAL_LOH) ||
               (somaticType() == Type.HET_DEL && tumorMinMinorAlleleCopies() < 0.5);
    }
}
