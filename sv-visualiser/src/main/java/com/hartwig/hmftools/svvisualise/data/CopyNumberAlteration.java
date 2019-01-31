package com.hartwig.hmftools.svvisualise.data;

import com.hartwig.hmftools.common.region.GenomeRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CopyNumberAlteration implements GenomeRegion {

    public abstract double copyNumber();

    public abstract double baf();

    public double minorAllelePloidy() {
        return Math.max(0, (1 - baf()) * copyNumber());
    }
}
