package com.hartwig.hmftools.datamodel.finding;

import java.util.List;

import com.hartwig.hmftools.datamodel.purple.PurpleMicrosatelliteStatus;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface MicrosatelliteStability extends Finding
{
    double microsatelliteIndelsPerMb();
    @NotNull PurpleMicrosatelliteStatus microsatelliteStatus();

    @NotNull List<LOHCopyNumbers> lohCopyNumbers();
}
