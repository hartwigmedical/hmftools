package com.hartwig.hmftools.serve.refgenome.liftover;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class LiftOverResult implements GenomePosition {

}
