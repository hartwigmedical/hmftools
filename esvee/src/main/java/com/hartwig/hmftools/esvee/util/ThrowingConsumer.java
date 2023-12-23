package com.hartwig.hmftools.esvee.util;

import java.util.function.Consumer;

public interface ThrowingConsumer<ARG> extends Consumer<ARG>
{
    static <ARG> Consumer<ARG> rethrow(final ThrowingConsumer<ARG> consumer)
    {
        return consumer;
    }

    void acceptThrowing(ARG argument) throws Exception;

    @Override
    default void accept(final ARG argument)
    {
        try
        {
            acceptThrowing(argument);
        }
        catch(final Exception e)
        {
            throw new RuntimeException(e);
        }
    }
}
