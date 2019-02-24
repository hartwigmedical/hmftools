package com.hartwig.hmftools.common.amber;

import java.util.EnumMap;
import java.util.Map;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface BaseDepth extends GenomePosition {

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

    int indelCount();

    default boolean isValid(int minAlleles) {
        return indelCount() == 0 && !ref().equals(Base.N) && baseMap().keySet().stream().filter(x -> !x.equals(Base.N)).count() >= minAlleles;
    }

    @NotNull
    default Base alt() {
        assert (isValid(1));

        Base result = Base.N;
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
        if (alt() == Base.N) {
            return 0;
        }

        return baseMap().getOrDefault(alt(), 0);
    }

}
