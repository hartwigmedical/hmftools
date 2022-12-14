package com.hartwig.hmftools.common.virus;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface VirusInterpreterData {

    @NotNull
    List<AnnotatedVirus> reportableViruses();

    // TODO Change to allViruses for consistency with other hmftools
    @NotNull
    List<AnnotatedVirus> unreportedViruses();
}
