package com.hartwig.hmftools.amber;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import com.hartwig.hmftools.common.amber.AmberBAF;

public class AmberUtils
{
    public static <T> List<T> getFuture(final List<Future<T>> futures) throws ExecutionException, InterruptedException
    {
        final List<T> result = new ArrayList<>();
        for(Future<T> chromosomeBAFEvidenceFuture : futures)
        {
            result.add(chromosomeBAFEvidenceFuture.get());
        }
        return result;
    }

    public static boolean isValid(final AmberBAF baf)
    {
        return Double.isFinite(baf.tumorBAF()) & Double.isFinite(baf.normalBAF());
    }
}
