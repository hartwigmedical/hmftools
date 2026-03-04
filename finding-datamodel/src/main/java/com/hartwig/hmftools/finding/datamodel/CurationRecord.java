package com.hartwig.hmftools.finding.datamodel;

import java.util.List;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record CurationRecord(
        @NotNull List<DriverCuration> driverCurations
)
{
    @SuppressWarnings("unused")
    public static JsonReadWriter<CurationRecord> jsonReadWriter()
    {
        return new JsonReadWriter<>(CurationRecord.class);
    }
}
