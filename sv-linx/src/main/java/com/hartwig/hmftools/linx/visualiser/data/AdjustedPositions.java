package com.hartwig.hmftools.linx.visualiser.data;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class AdjustedPositions
{

    public static List<AdjustedPosition> create(@NotNull final List<VisSvData> originalLinks, @NotNull final List<VisSvData> scaledLinks)
    {
        final List<AdjustedPosition> result = Lists.newArrayList();
        for (int i = 0; i < originalLinks.size(); i++)
        {
            final VisSvData original = originalLinks.get(i);
            final VisSvData scaled = scaledLinks.get(i);

            if (scaled.isValidStart())
            {
                final AdjustedPosition start = ImmutableAdjustedPosition.builder()
                        .chromosome(scaled.startChromosome())
                        .position(scaled.startPosition())
                        .unadjustedPosition(original.startPosition())
                        .svId(original.svId())
                        .build();

                result.add(start);
            }

            if (scaled.isValidEnd())
            {
                final AdjustedPosition end = ImmutableAdjustedPosition.builder()
                        .chromosome(scaled.endChromosome())
                        .position(scaled.endPosition())
                        .unadjustedPosition(original.endPosition())
                        .svId(original.svId())
                        .build();
                result.add(end);
            }
        }

        Collections.sort(result);
        return result;
    }
}
