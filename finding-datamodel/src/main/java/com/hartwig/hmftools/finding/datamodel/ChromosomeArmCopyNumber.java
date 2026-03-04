package com.hartwig.hmftools.finding.datamodel;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record ChromosomeArmCopyNumber(
        @NotNull String findingKey,
        @NotNull String chromosome,
        @NotNull ChromosomeArm arm,
        @NotNull Type type,
        double copyNumber,
        double relativeCopyNumber
) implements Finding
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
}
