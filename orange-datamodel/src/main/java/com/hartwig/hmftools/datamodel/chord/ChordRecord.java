package com.hartwig.hmftools.datamodel.chord;

import com.hartwig.hmftools.datamodel.finding.HomologousRecombination;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface ChordRecord
{
    // wrap the data to provide a finding type
    @NotNull
    HomologousRecombination homologousRecombination();

    // convenience accessors
    default double brca1Value()
    {
        return homologousRecombination().brca1Value();
    }
    default double brca2Value()
    {
        return homologousRecombination().brca2Value();
    }
    default double hrdValue()
    {
        return homologousRecombination().hrdValue();
    }
    @NotNull
    default ChordStatus hrStatus()
    {
        return homologousRecombination().hrStatus();
    }
    @NotNull
    default String hrdType()
    {
        return homologousRecombination().hrdType();
    }
}
