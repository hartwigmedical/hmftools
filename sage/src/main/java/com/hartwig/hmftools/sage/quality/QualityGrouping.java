package com.hartwig.hmftools.sage.quality;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class QualityGrouping {

    public static List<QualityCount> removePosition(final Collection<QualityCount> quality) {
        return groupBy(quality, QualityGrouping::removePosition);
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
        result.sort(QualityGrouping::compare);

        return result;
    }

    @NotNull
    private static QualityRecord removePosition(@NotNull final QualityCount count) {
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
    private static QualityRecord strandOnlyKey(@NotNull final QualityCount count) {
        return ImmutableQualityRecord.builder()
                .ref((byte) 'A')
                .alt((byte) 'A')
                .position(0)
                .qual((byte) 0)
                .trinucleotideContext()
                .readIndex(0)
                .firstOfPair(count.firstOfPair())
                .build();
    }

    private static int compare(final QualityCount o1, final QualityCount o2) {
        int countCompare = Integer.compare(o2.count(), o1.count());
        if (countCompare != 0) {
            return countCompare;
        }

        int refCompare = Byte.compare(o1.ref(), o2.ref());
        if (refCompare != 0) {
            return refCompare;
        }

        int altCompare = Byte.compare(o1.alt(), o2.alt());
        if (altCompare != 0) {
            return altCompare;
        }

        int indexCompare = Integer.compare(o1.readIndex(), o2.readIndex());
        if (indexCompare != 0) {
            return indexCompare;
        }

        int triOne = Byte.compare(o1.trinucleotideContext()[0], o2.trinucleotideContext()[0]);
        if (triOne != 0) {
            return triOne;
        }

        int triTwo = Byte.compare(o1.trinucleotideContext()[1], o2.trinucleotideContext()[1]);
        if (triTwo != 0) {
            return triTwo;
        }

        int triThree = Byte.compare(o1.trinucleotideContext()[2], o2.trinucleotideContext()[2]);
        if (triThree != 0) {
            return triThree;
        }

        return 0;
    }

}
