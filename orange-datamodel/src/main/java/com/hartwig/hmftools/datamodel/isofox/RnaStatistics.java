package com.hartwig.hmftools.datamodel.isofox;

import java.util.Set;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface RnaStatistics
{
    @NotNull
    Set<RnaQCStatus> qcStatus();

    long totalFragments();

    long duplicateFragments();

    double splicedFragmentPerc();
    double unsplicedFragmentPerc();
    double altFragmentPerc();
    double chimericFragmentPerc();
}
