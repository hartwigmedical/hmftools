package com.hartwig.hmftools.esvee.util;

import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Function;
import java.util.function.Predicate;

import org.jetbrains.annotations.Nullable;

public class Counter
{
    public final String Category;
    public final String Name;
    private final AtomicLong mValue = new AtomicLong();

    public Counter(final String name)
    {
        this("", name);
    }

    public Counter(final String category, final String name)
    {
        Category = category;
        Name = name;
    }

    public long getValue()
    {
        return mValue.get();
    }

    public void add(final long value)
    {
        mValue.addAndGet(value);
    }

    public <T, U> Function<T, U> wrap(final Function<T, U> function)
    {
        return argument ->
        {
            final long startTime = System.nanoTime();
            try
            {
                return function.apply(argument);
            }
            finally
            {
                final long duration = System.nanoTime() - startTime;
                add(duration);
            }
        };
    }

    public static <T, U> Function<T, U> wrap(@Nullable final Counter counter, final Function<T, U> function)
    {
        if(counter == null)
            return function;
        else
            return counter.wrap(function);
    }

    @Override
    public String toString()
    {
        return Name + ": " + getValue();
    }

    public String formatValue() { return String.valueOf(getValue()); }

    public static <T> Predicate<T> asPredicate(final Predicate<T> predicate, final Counter counter)
    {
        return item ->
        {
            if(predicate.test(item))
            {
                counter.add(1);
                return true;
            }
            return false;
        };
    }
}
