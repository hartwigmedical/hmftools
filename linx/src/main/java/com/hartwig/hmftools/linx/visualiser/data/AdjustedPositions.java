package com.hartwig.hmftools.linx.visualiser.data;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.visualiser.file.VisSvDataFile;

import org.jetbrains.annotations.NotNull;

public class AdjustedPositions
{

    public static List<AdjustedPosition> create(final List<VisSvDataFile> originalLinks, final List<VisSvDataFile> scaledLinks)
    {
        final List<AdjustedPosition> result = Lists.newArrayList();
        for (int i = 0; i < originalLinks.size(); i++)
        {
            final VisSvDataFile original = originalLinks.get(i);
            final VisSvDataFile scaled = scaledLinks.get(i);

            if (scaled.isValidStart())
            {
                final AdjustedPosition start = ImmutableAdjustedPosition.builder()
                        .chromosome(scaled.ChrStart)
                        .position(scaled.PosStart)
                        .unadjustedPosition(original.PosStart)
                        .svId(original.SvId)
                        .build();

                result.add(start);
            }

            if (scaled.isValidEnd())
            {
                final AdjustedPosition end = ImmutableAdjustedPosition.builder()
                        .chromosome(scaled.ChrEnd)
                        .position(scaled.PosEnd)
                        .unadjustedPosition(original.PosEnd)
                        .svId(original.SvId)
                        .build();
                result.add(end);
            }
        }

        Collections.sort(result);
        return result;
    }
}
