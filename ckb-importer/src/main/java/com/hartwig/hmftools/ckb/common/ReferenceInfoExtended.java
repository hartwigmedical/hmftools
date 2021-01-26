package com.hartwig.hmftools.ckb.common;

import org.jetbrains.annotations.Nullable;

public interface ReferenceInfoExtended {

    @Nullable
    public abstract String authors();

    @Nullable
    public abstract String journal();

    @Nullable
    public abstract String volume();

    @Nullable
    public abstract String issue();

    @Nullable
    public abstract String date();

    @Nullable
    public abstract String abstractText();

    @Nullable
    public abstract String year();
}
