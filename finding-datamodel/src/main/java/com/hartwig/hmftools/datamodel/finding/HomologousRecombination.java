package com.hartwig.hmftools.datamodel.finding;

import java.util.List;

import com.hartwig.hmftools.datamodel.chord.ChordStatus;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface HomologousRecombination extends Finding
{
    double brca1Value();
    double brca2Value();
    double hrdValue();

    @NotNull
    ChordStatus hrStatus();

    @NotNull
    String hrdType();

    @NotNull List<GainDeletion> lohCopyNumbers();

    @NotNull List<String> genes();
}
