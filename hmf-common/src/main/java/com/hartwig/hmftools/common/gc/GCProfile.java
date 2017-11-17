package com.hartwig.hmftools.common.gc;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface GCProfile extends GenomeRegion {

    double MIN_MAPPABLE_PERCENTAGE = 0.85;

    double gcContent();

    double nonNPercentage();

    double mappablePercentage();

    default boolean isMappable() {
        return Doubles.greaterOrEqual(mappablePercentage(), MIN_MAPPABLE_PERCENTAGE);
    }
}
