package com.hartwig.hmftools.common.purple.segment;

import java.util.List;

import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Cluster implements GenomeRegion {

    @NotNull
    public abstract List<GenomePosition> ratios();

    @NotNull
    public abstract List<StructuralVariantPosition> variants();

    @Nullable
    public Long firstVariant() {
        return variants().isEmpty() ? null : variants().get(0).position();
    }

    @Nullable
    public Long finalVariant() {
        return variants().isEmpty() ? null : variants().get(variants().size() - 1).position();
    }

    @Nullable
    public Long firstRatio() {
        return ratios().isEmpty() ? null : ratios().get(0).position();
    }

    @Nullable
    public Long finalRatio() {
        return ratios().isEmpty() ? null : ratios().get(ratios().size() - 1).position();
    }

    public StructuralVariantSupport type() {
        if (variants().isEmpty()) {
            return StructuralVariantSupport.NONE;
        }

        if (variants().size() == 1) {
            return StructuralVariantSupport.fromVariant(variants().get(0).type());
        }

        return StructuralVariantSupport.MULTIPLE;
    }
}
