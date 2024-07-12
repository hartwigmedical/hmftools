package com.hartwig.hmftools.common.sigs;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class SignatureTestFactory
{
    @NotNull
    public static ImmutableSignatureAllocation.Builder builder()
    {
        return ImmutableSignatureAllocation.builder().signature(Strings.EMPTY).allocation(0D).percent(0D);
    }
}
