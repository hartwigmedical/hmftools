package com.hartwig.hmftools.datamodel.virus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.List;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface VirusInterpreterData {

    @NotNull
    List<AnnotatedVirus> allViruses();

    @NotNull
    List<AnnotatedVirus> reportableViruses();
}
