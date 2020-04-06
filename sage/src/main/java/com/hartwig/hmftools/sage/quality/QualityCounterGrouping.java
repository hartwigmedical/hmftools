package com.hartwig.hmftools.sage.quality;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class QualityCounterGrouping {

    public static List<QualityCounter> groupByAlt(final Collection<QualityCounter> quality) {
        return groupBy(quality, QualityCounterGrouping::alt);
    }

    public static List<QualityCounter> groupWithoutPosition(final Collection<QualityCounter> quality) {
        return groupBy(quality, QualityCounterGrouping::withoutPosition);
    }

    @NotNull
    private static List<QualityCounter> groupBy(final Collection<QualityCounter> quality,
            final Function<QualityCounter, QualityCounterKey> keyFunction) {
        final Map<QualityCounterKey, QualityCounter> map = Maps.newHashMap();

        for (QualityCounter count : quality) {
            final QualityCounterKey key = keyFunction.apply(count);
            map.computeIfAbsent(key, QualityCounter::new).increment(count.count());
        }

        final List<QualityCounter> result = Lists.newArrayList(map.values());
        Collections.sort(result);

        return result;
    }

    @NotNull
    private static QualityCounterKey withoutPosition(@NotNull final QualityCounter count) {
        return ImmutableQualityCounterKey.builder().from(count).position(0).build();
    }

    @NotNull
    public static QualityCounterKey alt(@NotNull final QualityCounter count) {
        return ImmutableQualityCounterKey.builder().from(count).qual((byte) 0).trinucleotideContext().build();
    }

}
