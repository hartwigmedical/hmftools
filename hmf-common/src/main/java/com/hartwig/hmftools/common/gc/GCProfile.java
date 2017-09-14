package com.hartwig.hmftools.common.gc;

import com.hartwig.hmftools.common.region.GenomeRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface GCProfile extends GenomeRegion {

    double gcContent();

    double nonNPercentage();

    double mappablePercentage();
}
