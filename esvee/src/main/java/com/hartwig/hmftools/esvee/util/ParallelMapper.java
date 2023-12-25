package com.hartwig.hmftools.esvee.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Function;
import java.util.stream.Collectors;

import org.jetbrains.annotations.Nullable;

/** Why not Stream::parallel? Tasks that end up taking wildly differing time to process can leave quite a few idle threads
 * when using Stream::parallel. */
public final class ParallelMapper
{
    public static <ARG, RESULT> List<RESULT> flatMap(
            final ExecutorService executor, final List<ARG> items, final Function<ARG, List<RESULT>> mapper)
    {
        return map(executor, items, mapper).stream()
                .flatMap(Collection::stream)
                .collect(Collectors.toList());
    }

    public static <ARG, RESULT> List<RESULT> map(final ExecutorService executor, final List<ARG> items, final Function<ARG, RESULT> mapper)
    {
        return mapWithCounter(executor, items, mapper, null);
    }

    public static <ARG, RESULT> List<RESULT> mapWithCounter(
            final ExecutorService executor, final List<ARG> items, final Function<ARG, RESULT> mapper,
            @Nullable final Counter itemTimer)
    {
        final Function<ARG, RESULT> wrappedMapper = Counter.wrap(itemTimer, mapper);

        final List<CompletableFuture<RESULT>> futures = new ArrayList<>();
        for(final ARG ignored : items)
            futures.add(new CompletableFuture<>());

        final AtomicInteger index = new AtomicInteger();
        for (int i = 0; i < Runtime.getRuntime().availableProcessors(); i++)
        {
            executor.submit(() ->
            {
                while(true)
                {
                    final int toProcess = index.getAndIncrement();
                    if(toProcess >= items.size())
                        return;
                    final ARG item = items.get(toProcess);
                    final CompletableFuture<RESULT> future = futures.get(toProcess);
                    try
                    {
                        future.complete(wrappedMapper.apply(item));
                    }
                    catch(final Throwable throwable)
                    {
                        future.completeExceptionally(throwable);
                    }
                }
            });
        }

        return futures.stream()
                .map(ThrowingFunction.rethrow(Future::get))
                .filter(Objects::nonNull)
                .collect(Collectors.toList());
    }

    public static <ARG, RESULT> List<RESULT> mapWithProgress(
            final ExecutorService executor, final List<ARG> items, final Function<ARG, RESULT> mapper)
    {
            return map(executor, items, item -> mapper.apply(item));
    }
}
