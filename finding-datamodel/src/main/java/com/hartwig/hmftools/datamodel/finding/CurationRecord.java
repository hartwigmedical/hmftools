package com.hartwig.hmftools.datamodel.finding;

import java.util.List;

import io.soabase.recordbuilder.core.RecordBuilder;
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
