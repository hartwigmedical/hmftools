package com.hartwig.hmftools.geneutils.probequality;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Stream;

public class Utils
{
    // Splits the stream into lists of at most `partitionSize` elements.
    public static <T> Stream<List<T>> partitionStream(Stream<T> elements, int partitionSize)
    {
        if(partitionSize < 1)
        {
            throw new IllegalArgumentException("partitionSize must be >= 1");
        }

        Iterator<T> iterator = elements.iterator();
        return Stream.generate(() ->
        {
            List<T> batch = new ArrayList<>(partitionSize);
            for(int i = 0; i < partitionSize && iterator.hasNext(); i++)
            {
                batch.add(iterator.next());
            }
            return batch;
        }).takeWhile(b -> !b.isEmpty());
    }
}
