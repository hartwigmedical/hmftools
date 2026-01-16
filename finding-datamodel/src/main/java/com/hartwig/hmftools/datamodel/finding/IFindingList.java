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
    @NotNull List<T> all();

    default int size()
    {
        return all().size();
    }

    default boolean isEmpty()
    {
        return all().isEmpty();
    }

    default Iterator<T> iterator()
    {
        return all().iterator();
    }

    default T get(int index)
    {
        return all().get(index);
    }

    default ListIterator<T> listIterator()
    {
        return all().listIterator();
    }

    default ListIterator<T> listIterator(int index)
    {
        return all().listIterator(index);
    }

    default Stream<T> stream()
    {
        return all().stream();
    }

    default Stream<T> parallelStream()
    {
        return all().parallelStream();
    }
}
