package com.hartwig.hmftools.sage.quality;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class QualityGrouping {

    public static Collection<QualityCount> groupWithoutPosition(final Collection<QualityCount> quality) {
        return groupBy(quality, QualityGrouping::withoutPosition);
    }

    public static Collection<QualityCount> groupByQuality(final Collection<QualityCount> quality) {
        return groupBy(quality, QualityGrouping::qualAccuracyKey);
    }

    public static Collection<QualityCount> groupByStrand(final Collection<QualityCount> quality) {
        return groupBy(quality, QualityGrouping::strandKey);
    }

    public static Collection<QualityCount> groupByStrandOnly(final Collection<QualityCount> quality) {
        return groupBy(quality, QualityGrouping::strandOnlyKey);
    }

    @NotNull
    private static Collection<QualityCount> groupBy(final Collection<QualityCount> quality,
            final Function<QualityCount, QualityRecord> keyFunction) {
        final Map<QualityRecord, QualityCount> map = Maps.newHashMap();

        for (QualityCount count : quality) {
            final QualityRecord key = keyFunction.apply(count);
            map.computeIfAbsent(key, QualityCount::new).increment(count.count());
        }

        final List<QualityCount> result = Lists.newArrayList(map.values());
        result.sort(Comparator.comparingInt(QualityCount::count).reversed());

        return result;
    }

    @NotNull
    private static QualityRecord withoutPosition(@NotNull final QualityCount count) {
        return ImmutableQualityRecord.builder()
                .ref(count.ref())
                .alt(count.alt())
                .qual(count.qual())
                .firstOfPair(count.firstOfPair())
                .position(0)
                .build();
    }

    @NotNull
    private static QualityRecord qualAccuracyKey(@NotNull final QualityCount count) {
        return ImmutableQualityRecord.builder().ref(count.ref()).alt(count.alt()).qual(count.qual()).firstOfPair(false).position(0).build();
    }

    @NotNull
    private static QualityRecord strandKey(@NotNull final QualityCount count) {
        return ImmutableQualityRecord.builder()
                .ref(count.ref())
                .alt(count.alt())
                .position(0)
                .qual((byte) 0)
                .firstOfPair(count.firstOfPair())
                .build();
    }

    @NotNull
    private static QualityRecord strandOnlyKey(@NotNull final QualityCount count) {
        return ImmutableQualityRecord.builder()
                .ref((byte) 'A')
                .alt((byte) 'A')
                .position(0)
                .qual((byte) 0)
                .firstOfPair(count.firstOfPair())
                .build();
    }

}
