package com.hartwig.hmftools.sage.quality;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;

public class QualityRecalibrationMap
{

    @NotNull
    private final Map<QualityRecalibrationKey, QualityRecalibrationRecord> map;

    public QualityRecalibrationMap(@NotNull final List<QualityRecalibrationRecord> records)
    {
        this.map = records.stream().collect(Collectors.toMap(QualityRecalibrationRecord::key, x -> x));
    }

    public double quality(byte ref, byte alt, byte[] trinucleotideContext, byte qual)
    {
        final QualityRecalibrationKey key =
                ImmutableQualityRecalibrationKey.builder().ref(ref).alt(alt).qual(qual).trinucleotideContext(trinucleotideContext).build();

        return Optional.ofNullable(map.get(key)).map(QualityRecalibrationRecord::recalibratedQual).orElse(qual * 1d);
    }
}
