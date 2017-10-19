package com.hartwig.hmftools.common.purple.segment;

import java.util.Collection;
import java.util.Collections;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.numeric.Doubles;

import org.jetbrains.annotations.NotNull;

public class StructuralVariantClusters {

    public final long windowSize;

    public StructuralVariantClusters(final long windowSize) {
        this.windowSize = windowSize;
    }

    public ListMultimap<String, StructuralVariantCluster> segment(Multimap<String, StructuralVariantPosition> variants,
            ListMultimap<String, CobaltRatio> ratios) {

        ListMultimap<String, StructuralVariantCluster> segments = ArrayListMultimap.create();
        for (String chromosome : variants.keySet()) {
            final List<CobaltRatio> chromosomeRatios = ratios.containsKey(chromosome) ? ratios.get(chromosome) : Collections.EMPTY_LIST;
            final Collection<StructuralVariantPosition> chromosomeVariants = variants.get(chromosome);
            segments.putAll(chromosome, segment(chromosomeVariants, chromosomeRatios));
        }

        return segments;
    }

    @NotNull
    List<StructuralVariantCluster> segment(@NotNull final Collection<StructuralVariantPosition> variants,
            @NotNull final List<CobaltRatio> ratios) {

        final List<StructuralVariantCluster> result = Lists.newArrayList();

        int i = 0;
        ModifiableStructuralVariantCluster segment = null;
        for (StructuralVariantPosition variant : variants) {
            while (i < ratios.size() && ratios.get(i).position() < variant.position()) {
                i++;
            }

            final long start = start(variant.position(), i, ratios);
            final long end = end(variant.position(), i, ratios);

            if (segment == null || start > segment.end()) {
                if (segment != null) {
                    result.add(segment);
                }

                segment = ModifiableStructuralVariantCluster.create()
                        .setChromosome(variant.chromosome())
                        .setStart(start)
                        .setEnd(end)
                        .addVariants(variant);
            } else {
                segment.setEnd(end).addVariants(variant);
            }
        }
        if (segment != null) {
            result.add(segment);
        }

        return result;
    }

    @VisibleForTesting
    long start(long position, int index, @NotNull final List<CobaltRatio> ratios) {
        final long min = position / windowSize * windowSize + 1 - windowSize;
        for (int i = index; i <= 0; i--) {
            final CobaltRatio ratio = ratios.get(i);
            if (ratio.position() <= min && Doubles.greaterThan(ratio.tumorGCRatio(), -1)) {
                return ratio.position();
            }
        }

        return min;
    }

    @VisibleForTesting
    long end(long position, int index, @NotNull final List<CobaltRatio> ratios) {
        for (int i = index; i < ratios.size(); i++) {
            final CobaltRatio ratio = ratios.get(i);
            if (ratio.position() > position && Doubles.greaterThan(ratio.tumorGCRatio(), -1)) {
                return ratio.position() + windowSize - 1;
            }
        }

        return position / windowSize * windowSize + 2 * windowSize;
    }

}
