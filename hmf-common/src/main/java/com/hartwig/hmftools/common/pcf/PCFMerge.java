package com.hartwig.hmftools.common.pcf;

import java.util.Comparator;
import java.util.List;
import java.util.function.Function;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

public class PCFMerge {

    private static final Comparator<PCFPosition> COMPARE =
            Comparator.comparing((Function<PCFPosition, String>) GenomePosition::chromosome).thenComparingLong(GenomePosition::position);

    @NotNull
    static List<ModifiablePCFPosition> merge(@NotNull final List<ModifiablePCFPosition> positions) {

        positions.sort(COMPARE);

        int i = 0;
        while (i < positions.size() - 1) {
            ModifiablePCFPosition current = positions.get(i);
            ModifiablePCFPosition next = positions.get(i + 1);

            if (current.chromosome().equals(next.chromosome()) && current.position() == next.position()) {
                current.setMinPosition(Math.max(current.minPosition(), next.minPosition()))
                        .setMaxPosition(Math.min(current.maxPosition(), next.maxPosition()));
                positions.remove(i + 1);
            } else if (current.chromosome().equals(next.chromosome())) {
                current.setMaxPosition(Math.min(current.maxPosition(), next.position()));
                next.setMinPosition(Math.max(next.minPosition(), current.position()));
                i++;
            } else {
                i++;
            }
        }

        return positions;
    }

}
