package com.hartwig.hmftools.svassembly.models;

import org.jetbrains.annotations.Nullable;

public interface TrimmableAssembly<SELF>
{
    @Nullable
    SELF trim(final int removeLeft, final int removeRight);
}
