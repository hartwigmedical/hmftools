package com.hartwig.hmftools.esvee.util;

import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.Supplier;

import org.jetbrains.annotations.Nullable;

public class Counter
{
    public final String Category;
    public final String Name;
    public final boolean IsTime;
    private final AtomicLong mValue = new AtomicLong();

    public Counter(final String name)
    {
        this(name, false);
    }

    public Counter(final String name, final boolean isTime)
    {
        this("", name, isTime);
    }

    public Counter(final String category, final String name)
    {
        this(category, name, false);
    }

    public Counter(final String category, final String name, final boolean isTime)
    {
        Category = category;
        Name = name;
        IsTime = isTime;
    }

    public long getValue()
    {
        return mValue.get();
    }

    public void add(final long value)
    {
        mValue.addAndGet(value);
    }

    public void time(final Runnable runnable)
    {
        final long startTime = System.nanoTime();
        try
        {
            runnable.run();
        }
        finally
        {
            final long duration = System.nanoTime() - startTime;
            add(duration);
        }
    }

    public <T> T time(final Supplier<T> supplier)
    {
        final long startTime = System.nanoTime();
        try
        {
            return supplier.get();
        }
        finally
        {
            final long duration = System.nanoTime() - startTime;
            add(duration);
        }
    }

    public static <T> T time(@Nullable final Counter counter, final Supplier<T> supplier)
    {
        if(counter == null)
            return supplier.get();
        else
            return counter.time(supplier);
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

    /*
    previous implemented interface CSVWriter.CSVValue asValueForCSV
     */

    public Object asValueForCSV()
    {
        return getValue();
    }

    public String formatValue()
    {
        if (IsTime)
            return CommonUtils.formatNanos(getValue());
        else
            return String.valueOf(getValue());
    }

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
