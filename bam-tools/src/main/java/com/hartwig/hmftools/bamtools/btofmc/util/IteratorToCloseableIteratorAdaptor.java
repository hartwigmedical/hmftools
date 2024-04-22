package com.hartwig.hmftools.bamtools.btofmc.util;

import java.util.Iterator;

import htsjdk.samtools.util.CloseableIterator;

public class IteratorToCloseableIteratorAdaptor<T> implements CloseableIterator<T>
{
    private final Iterator<T> mIterator;

    public IteratorToCloseableIteratorAdaptor(Iterator<T> iterator)
    {
        mIterator = iterator;
    }

    @Override
    public void close()
    {
    }

    @Override
    public boolean hasNext()
    {
        return mIterator.hasNext();
    }

    @Override
    public T next()
    {
        return mIterator.next();
    }
}
