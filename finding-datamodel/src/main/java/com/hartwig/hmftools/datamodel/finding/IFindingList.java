package com.hartwig.hmftools.datamodel.finding;

import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.stream.Stream;

import jakarta.validation.constraints.NotNull;

@SuppressWarnings("unused")
public interface IFindingList<T extends Finding>
{
    @NotNull FindingsStatus status();
    @NotNull List<T> findings();

    default int size()
    {
        return findings().size();
    }

    default boolean isEmpty()
    {
        return findings().isEmpty();
    }

    default Iterator<T> iterator()
    {
        return findings().iterator();
    }

    default T get(int index)
    {
        return findings().get(index);
    }

    default ListIterator<T> listIterator()
    {
        return findings().listIterator();
    }

    default ListIterator<T> listIterator(int index)
    {
        return findings().listIterator(index);
    }

    default Stream<T> stream()
    {
        return findings().stream();
    }

    default Stream<T> parallelStream()
    {
        return findings().parallelStream();
    }
}
