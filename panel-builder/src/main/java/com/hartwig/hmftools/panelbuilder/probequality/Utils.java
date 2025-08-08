package com.hartwig.hmftools.panelbuilder.probequality;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Stream;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;

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
            ArrayList<T> batch = new ArrayList<>(partitionSize);
            for(int i = 0; i < partitionSize && iterator.hasNext(); i++)
            {
                batch.add(iterator.next());
            }
            return (List<T>) batch;
        }).takeWhile(b -> !b.isEmpty());
    }

    public static BwaMemAligner createBwaMemAligner(String bwaIndexImageFile, int threads)
    {
        if(!Files.exists(Paths.get(bwaIndexImageFile)) || bwaIndexImageFile.isEmpty())
        {
            throw new RuntimeException("Reference genome file is missing or empty");
        }

        BwaMemIndex index = new BwaMemIndex(bwaIndexImageFile);
        BwaMemAligner aligner = new BwaMemAligner(index);

        aligner.setNThreadsOption(threads);

        return aligner;
    }
}
