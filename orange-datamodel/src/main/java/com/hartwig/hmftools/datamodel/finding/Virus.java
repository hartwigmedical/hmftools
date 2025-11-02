package com.hartwig.hmftools.datamodel.finding;

import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface Virus extends Finding, VirusInterpreterEntry {
}
