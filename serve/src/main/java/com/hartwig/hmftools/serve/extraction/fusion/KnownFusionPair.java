package com.hartwig.hmftools.serve.extraction.fusion;

import com.hartwig.hmftools.common.serve.datamodel.fusion.FusionPair;
import com.hartwig.hmftools.serve.extraction.KnownEvent;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class KnownFusionPair implements FusionPair, KnownEvent {

}
