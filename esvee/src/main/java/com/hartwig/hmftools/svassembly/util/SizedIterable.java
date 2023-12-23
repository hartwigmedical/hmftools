package com.hartwig.hmftools.svassembly.util;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Spliterators;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import org.jetbrains.annotations.NotNull;

public class SizedIterable<T> implements Iterable<T>
{
    private int mSize;
    private final Supplier<Iterator<T>> mIteratorFactory;

    public SizedIterable(final int size, final Supplier<Iterator<T>> iteratorFactory)
    {
        mSize = size;
        mIteratorFactory = iteratorFactory;
    }

    public static <T, U> SizedIterable<T> create(final int size, final Supplier<Iterator<U>> iteratorFactory, final Function<U, T> converter)
    {
        return new SizedIterable<>(size, new Supplier<Iterator<T>>()
        {
            @Override
            public Iterator<T> get()
            {
                return new Iterator<>()
                {
                    final Iterator<U> inner = iteratorFactory.get();

                    @Override
                    public boolean hasNext()
                    {
                        return inner.hasNext();
                    }

                    @Override
                    public T next()
                    {
                        return converter.apply(inner.next());
                    }
                };
            }
        });
    }

    public int size()
    {
        if (mSize == -1)
            mSize = (int) stream().count();
        return mSize;
    }

    @NotNull
    @Override
    public Iterator<T> iterator()
    {
        return mIteratorFactory.get();
    }

    public Stream<T> stream()
    {
        return StreamSupport.stream(Spliterators.spliterator(iterator(), mSize, 0), false);
    }

    public List<T> toList()
    {
        final List<T> list = new ArrayList<>(mSize);
        for (final T item : this)
            list.add(item);
        return list;
    }
}
