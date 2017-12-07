package com.hartwig.hmftools.common.purple.segment;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.pcf.PCFPosition;
import com.hartwig.hmftools.common.pcf.PCFSource;
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
    public abstract List<PCFPosition> pcfPositions();

    @NotNull
    public abstract List<ClusterVariantLeg> variants();

    @NotNull
    public List<GenomePosition> ratios() {
        return pcfPositions().stream().filter(x -> !x.source().equals(PCFSource.TUMOR_BAF)).collect(Collectors.toList());
    }


    @Nullable
    Long firstVariant() {
        return variants().isEmpty() ? null : variants().get(0).position();
    }

    @Nullable
    Long finalVariant() {
        return variants().isEmpty() ? null : variants().get(variants().size() - 1).position();
    }

    @Nullable
    Long firstRatio() {
        return ratios().isEmpty() ? null : ratios().get(0).position();
    }

    @Nullable
    Long finalRatio() {
        return ratios().isEmpty() ? null : ratios().get(ratios().size() - 1).position();
    }

}
