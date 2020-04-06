package com.hartwig.hmftools.sage.quality;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class QualityRecalibrationFactory {

    @NotNull
    public static List<QualityRecalibrationRecord> create(@NotNull final List<QualityCounter> records) {

        final List<QualityRecalibrationRecord> result = Lists.newArrayList();
        final List<QualityCounter> cleanedRecords = QualityCounterGrouping.groupWithoutPosition(records)
                .stream()
                .filter(x -> (char) x.ref() != 'N')
                .filter(x -> (char) x.alt() != 'N')
                .collect(Collectors.toList());

        final Map<QualityRecalibrationKey, Integer> refCountMap = cleanedRecords.stream()
                .filter(x -> x.ref() == x.alt())
                .collect(Collectors.toMap(QualityRecalibrationFactory::key, QualityCounter::count));

        for (QualityCounter cleanedRecord : cleanedRecords) {
            final QualityRecalibrationKey key = key(cleanedRecord);
            final QualityRecalibrationKey refKey = ImmutableQualityRecalibrationKey.builder().from(key).alt(cleanedRecord.ref()).build();
            int refCount = refCountMap.getOrDefault(refKey, 0);
            if (refCount > 0) {
                double recalibratedQual = cleanedRecord.alt() == cleanedRecord.ref()
                        ? cleanedRecord.qual()
                        : recalibratedQual(refCount, cleanedRecord.count());

                result.add(ImmutableQualityRecalibrationRecord.builder()
                        .key(key)
                        .count(cleanedRecord.count())
                        .recalibratedQual(recalibratedQual)
                        .build());

            }
        }

        return result;
    }

    @VisibleForTesting
    static double recalibratedQual(int refCount, int altCount) {
        double percent = altCount / (double) refCount;
        return -10 * Math.log10(percent);
    }

    @NotNull
    static QualityRecalibrationKey key(@NotNull final QualityCounterKey record) {
        return ImmutableQualityRecalibrationKey.builder()
                .ref(record.ref())
                .alt(record.alt())
                .qual(record.qual())
                .trinucleotideContext(record.trinucleotideContext())
                .build();
    }

}
