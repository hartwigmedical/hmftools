package com.hartwig.hmftools.datamodel.purple;

import com.hartwig.hmftools.datamodel.finding.MicrosatelliteStability;
import com.hartwig.hmftools.datamodel.finding.TumorMutationStatus;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PurpleCharacteristics
{
    boolean wholeGenomeDuplication();

    @NotNull
    MicrosatelliteStability microsatelliteStability();

    @NotNull
    TumorMutationStatus tumorMutationStatus();
}
