package com.hartwig.hmftools.common.variant;

import com.hartwig.hmftools.common.position.GenomePosition;
import org.jetbrains.annotations.NotNull;

public interface Variant extends GenomePosition {

    @NotNull
    String ref();

    @NotNull
    String alt();

    @NotNull
    VariantType type();

    @NotNull
    String filter();
}
