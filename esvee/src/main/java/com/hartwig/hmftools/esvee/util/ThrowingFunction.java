package com.hartwig.hmftools.esvee.util;

import java.util.function.Function;

public interface ThrowingFunction<ARG, RET> extends Function<ARG, RET>
{
    static <ARG, RET> Function<ARG, RET> rethrow(final ThrowingFunction<ARG, RET> function)
    {
        return function;
    }

    RET applyThrowing(ARG argument) throws Exception;

    @Override
    default RET apply(final ARG argument)
    {
        try
        {
            return applyThrowing(argument);
        }
        catch(final Exception e)
        {
            throw new RuntimeException(e);
        }
    }
}
