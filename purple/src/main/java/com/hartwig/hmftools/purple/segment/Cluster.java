package com.hartwig.hmftools.purple.segment;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.common.utils.pcf.PCFSource;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
abstract class Cluster implements GenomeRegion
{
    @NotNull
    public abstract List<PCFPosition> pcfPositions();

    @NotNull
    public abstract List<SVSegment> variants();

    @NotNull
    public List<GenomePosition> ratios()
    {
        return pcfPositions().stream().filter(x -> !x.source().equals(PCFSource.TUMOR_BAF)).collect(Collectors.toList());
    }

}
