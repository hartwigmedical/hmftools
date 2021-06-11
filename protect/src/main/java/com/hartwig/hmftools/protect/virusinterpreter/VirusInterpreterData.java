package com.hartwig.hmftools.protect.virusinterpreter;

import java.util.List;

import com.hartwig.hmftools.common.virus.AnnotatedVirus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface VirusInterpreterData {

    @NotNull
    List<AnnotatedVirus> reportableViruses();

    @NotNull
    List<AnnotatedVirus> unreportedViruses();
}
