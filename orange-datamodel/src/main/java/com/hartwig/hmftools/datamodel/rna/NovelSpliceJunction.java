package com.hartwig.hmftools.datamodel.rna;

import org.immutables.value.Value;

@Value.Immutable
public interface NovelSpliceJunction {
    String geneName();

    String chromosome();
    int junctionStart();
    int junctionEnd();
    AltSpliceJunctionType type();
    int fragmentCount();
    int depthStart();
    int depthEnd();
    AltSpliceJunctionContext regionStart();
    AltSpliceJunctionContext regionEnd();
    int cohortFrequency();
}
