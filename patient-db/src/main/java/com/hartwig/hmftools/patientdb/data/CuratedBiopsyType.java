package com.hartwig.hmftools.patientdb.data;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class CuratedBiopsyType {

    @Nullable
    public abstract String type();

    @Nullable
    public abstract String searchPrimaryTumorLocation();

    @Nullable
    public abstract String searchCancerSubType();

    @Nullable
    public abstract String searchBiopsySite();

    @Nullable
    public abstract String searchBiopsyLocation();
}
