package com.hartwig.hmftools.esvee.sequence;

import org.jetbrains.annotations.Nullable;

public interface TrimmableAssembly<SELF>
{
    @Nullable
    SELF trim(final int removeLeft, final int removeRight);
}
