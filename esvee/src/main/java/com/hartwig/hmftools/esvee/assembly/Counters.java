package com.hartwig.hmftools.esvee.assembly;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Streams;
import com.hartwig.hmftools.esvee.util.Counter;
import com.hartwig.hmftools.esvee.util.ThrowingFunction;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public abstract class Counters<T extends Counters<T>>
{
    private static final Logger LOGGER = LogManager.getLogger(Counters.class);
    private static final Map<Class<?>, Function<Object, List<Counter>>> COUNTER_ACCESSORS = new ConcurrentHashMap<>();

    @SuppressWarnings("unused")
    public void log()
    {
        if (LOGGER.isTraceEnabled())
            all().forEach(counter -> LOGGER.trace("{}", counter));
    }

    public List<Counter> all()
    {
        return COUNTER_ACCESSORS.computeIfAbsent(getClass(), Counters::createAccessor).apply(this);
    }

    public Counters<T> add(final Counters<T> other)
    {
        if (other.getClass() != getClass())
            throw new IllegalArgumentException("Cannot add non-equal counter types "
                    + getClass().getSimpleName() + " and " + other.getClass().getSimpleName());

        // We guarantee a deterministic iteration order, so we can just zip these.
        //noinspection UnstableApiUsage
        Streams.forEachPair(all().stream(), other.all().stream(), (l, r) -> l.add(r.getValue()));
        return this;
    }

    private static Function<Object, List<Counter>> createAccessor(final Class<?> clazz)
    {
        final List<ThrowingFunction<Object, Counter>> accessors = new ArrayList<>();
        for(final Field field : clazz.getDeclaredFields())
        {
            if(!field.getType().equals(Counter.class))
                continue;

            field.setAccessible(true);
            accessors.add(instance -> (Counter) (field.get(instance)));
        }

        return instance -> accessors.stream()
                .map(accessor -> accessor.apply(instance))
                .collect(Collectors.toList());
    }
}
