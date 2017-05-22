package com.hartwig.hmftools.common.region.bed;

import com.hartwig.hmftools.common.region.GenomeRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true, passAnnotations = { NotNull.class, Nullable.class})
public abstract class BEDGenomeRegion implements GenomeRegion {
}
