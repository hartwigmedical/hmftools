package com.hartwig.hmftools.sage.context;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface RefContextCandidates {

    @Nullable
    RefContext refContext(@NotNull final String chromosome, final long position);
}
