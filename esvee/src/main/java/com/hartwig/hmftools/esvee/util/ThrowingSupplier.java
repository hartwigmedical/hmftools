package com.hartwig.hmftools.esvee.util;

import java.util.function.Supplier;

public interface ThrowingSupplier<RET> extends Supplier<RET>
{
    static <RET> Supplier<RET> rethrow(final ThrowingSupplier<RET> supplier)
    {
        return supplier;
    }

    RET getThrowing() throws Exception;

    @Override
    default RET get()
    {
        try
        {
            return getThrowing();
        }
        catch(final Exception e)
        {
            throw new RuntimeException(e);
        }
    }
}
