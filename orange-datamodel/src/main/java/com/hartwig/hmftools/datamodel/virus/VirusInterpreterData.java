package com.hartwig.hmftools.datamodel.virus;

import java.util.List;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface VirusInterpreterData
{
    @NotNull
    List<VirusInterpreterEntry> allViruses();

    @Gson.Ignore
    @NotNull
    default List<VirusInterpreterEntry> reportableViruses()
    {
        return allViruses().stream().filter(VirusInterpreterEntry::reported).toList();
    }
}
