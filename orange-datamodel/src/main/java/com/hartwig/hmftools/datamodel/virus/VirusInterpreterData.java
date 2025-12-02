package com.hartwig.hmftools.datamodel.virus;

import java.util.List;

import com.hartwig.hmftools.datamodel.finding.Virus;

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

    @NotNull
    List<Virus> driverViruses();
}
