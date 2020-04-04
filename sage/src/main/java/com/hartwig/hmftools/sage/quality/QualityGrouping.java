package com.hartwig.hmftools.sage.quality;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class QualityGrouping {

    public static List<QualityCount> groupByAlt(final Collection<QualityCount> quality) {
        return groupBy(quality, QualityGrouping::alt);
    }

    public static List<QualityCount> groupWithoutPosition(final Collection<QualityCount> quality) {
        return groupBy(quality, QualityGrouping::withoutPosition);
    }

    public static List<QualityCount> groupByQuality(final Collection<QualityCount> quality) {
        return groupBy(quality, QualityGrouping::qualAccuracyKey);
    }

    public static List<QualityCount> groupByStrand(final Collection<QualityCount> quality) {
        return groupBy(quality, QualityGrouping::strandKey);
    }

    public static List<QualityCount> groupByStrandOnly(final Collection<QualityCount> quality) {
        return groupBy(quality, QualityGrouping::strandOnlyKey);
    }

    @NotNull
    private static List<QualityCount> groupBy(final Collection<QualityCount> quality,
            final Function<QualityCount, QualityRecord> keyFunction) {
        final Map<QualityRecord, QualityCount> map = Maps.newHashMap();

        for (QualityCount count : quality) {
            final QualityRecord key = keyFunction.apply(count);
            map.computeIfAbsent(key, QualityCount::new).increment(count.count());
        }

        final List<QualityCount> result = Lists.newArrayList(map.values());
        Collections.sort(result);

        return result;
    }

    @NotNull
    private static QualityRecord withoutPosition(@NotNull final QualityCount count) {
        return ImmutableQualityRecord.builder().from(count).position(0).build();
    }

    @NotNull
    private static QualityRecord qualAccuracyKey(@NotNull final QualityCount count) {
        return ImmutableQualityRecord.builder().from(count).position(0).firstOfPair(false).build();
    }

    @NotNull
    private static QualityRecord strandKey(@NotNull final QualityCount count) {
        return ImmutableQualityRecord.builder().from(count).position(0).qual((byte) 0).build();
    }

    @NotNull
    public static QualityRecord alt(@NotNull final QualityCount count) {
        return ImmutableQualityRecord.builder().from(count).qual((byte) 0).trinucleotideContext().firstOfPair(false).build();
    }

    @NotNull
    private static QualityRecord strandOnlyKey(@NotNull final QualityCount count) {
        return ImmutableQualityRecord.builder()
                .ref((byte) 'A')
                .alt((byte) 'A')
                .position(0)
                .qual((byte) 0)
                .trinucleotideContext()
                .firstOfPair(count.firstOfPair())
                .build();
    }

}
