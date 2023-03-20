package com.hartwig.hmftools.datamodel.rna;

import org.immutables.value.Value;

@Value.Immutable
public interface RnaFusion {
    String name();
    String chromosomeUp();
    String chromosomeDown();
    int positionUp();
    int positionDown();
    String junctionTypeUp();
    String junctionTypeDown();
    StructuralVariantType svType();
    int splitFragments();
    int realignedFrags();
    int discordantFrags();
    int depthUp();
    int depthDown();
    int cohortFrequency();
}
