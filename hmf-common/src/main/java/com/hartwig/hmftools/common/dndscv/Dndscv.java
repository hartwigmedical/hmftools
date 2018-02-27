package com.hartwig.hmftools.common.dndscv;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface Dndscv {

    @NotNull
    String gene();

    int nSynonymous();

    int nMissense();

    int nNonSynonymous();

    int nSplice();

    int nIndel();

    double wMissense();

    double wNonSynonymous();

    double wSplice();

    double wIndel();

    double pScore();

    double qScore();
}
