package com.hartwig.hmftools.common.amber;

import java.util.EnumMap;
import java.util.Map;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface NormalBAF extends GenomePosition {

    enum Base {
        G,
        A,
        T,
        C,
        N
    }

    @NotNull
    Base ref();

    @NotNull
    EnumMap<Base, Integer> baseMap();

    int readDepth();

    default boolean isValid() {
        return !ref().equals(Base.N) && baseMap().keySet().stream().filter(x -> !x.equals(Base.N)).count() > 1;
    }

    @NotNull
    default Base alt() {
        assert (isValid());

        Base result = ref();
        int maxCount = 0;
        for (Map.Entry<Base, Integer> entry : baseMap().entrySet()) {
            if (!entry.getKey().equals(ref()) && entry.getValue() > maxCount) {
                maxCount = entry.getValue();
                result = entry.getKey();
            }

        }
        return result;
    }

    default int refSupport() {
        return baseMap().getOrDefault(ref(), 0);
    }

    default int altSupport() {
        return baseMap().getOrDefault(alt(), 0);
    }

}
