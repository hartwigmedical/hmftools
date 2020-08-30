package com.hartwig.hmftools.serve.hartwigcohort;

import org.jetbrains.annotations.NotNull;

public interface HartwigEntry {

    @NotNull
    String chromosome();

    long position();

    @NotNull
    String ref();

    @NotNull
    String alt();

    @NotNull
    String gene();

    @NotNull
    String transcript();

    @NotNull
    String proteinAnnotation();
}
