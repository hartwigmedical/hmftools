package com.hartwig.hmftools.common.rna;

import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.immutables.value.Value;

@Value.Immutable
public abstract class RnaFusion
{
    public abstract String name();
    public abstract String chromosomeUp();
    public abstract String chromosomeDown();
    public abstract int positionUp();
    public abstract int positionDown();
    public abstract byte orientationUp();
    public abstract byte orientationDown();
    public abstract String junctionTypeUp();
    public abstract String junctionTypeDown();
    public abstract StructuralVariantType svType();
    public abstract int splitFragments();
    public abstract int realignedFrags();
    public abstract int discordantFrags();
    public abstract int depthUp();
    public abstract int depthDown();
    public abstract int maxAnchorLengthUp();
    public abstract int maxAnchorLengthDown();
    public abstract int cohortFrequency();
}
