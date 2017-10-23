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
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.jetbrains.annotations.NotNull;

public class StructuralVariantClusterFactory {

    public final long windowSize;

    public StructuralVariantClusterFactory(final long windowSize) {
        this.windowSize = windowSize;
    }

    @NotNull
    public ListMultimap<String, StructuralVariantCluster> cluster(@NotNull final List<StructuralVariant> variants,
            ListMultimap<String, CobaltRatio> ratios) {
        final Multimap<String, StructuralVariantPosition> positions = asMap(StructuralVariantPositionFactory.create(variants));
        return cluster(positions, ratios);
    }

    @NotNull
    public ListMultimap<String, StructuralVariantCluster> cluster(@NotNull final Multimap<String, StructuralVariantPosition> variants,
            @NotNull final ListMultimap<String, CobaltRatio> ratios) {

        ListMultimap<String, StructuralVariantCluster> segments = ArrayListMultimap.create();
        for (String chromosome : variants.keySet()) {
            final List<CobaltRatio> chromosomeRatios = ratios.containsKey(chromosome) ? ratios.get(chromosome) : Collections.EMPTY_LIST;
            final Collection<StructuralVariantPosition> chromosomeVariants = variants.get(chromosome);
            segments.putAll(chromosome, cluster(chromosomeVariants, chromosomeRatios));
        }

        return segments;
    }

    @NotNull
    List<StructuralVariantCluster> cluster(@NotNull final Collection<StructuralVariantPosition> variants,
            @NotNull final List<CobaltRatio> ratios) {

        final List<StructuralVariantCluster> result = Lists.newArrayList();

        int ratioIndex = 0;
        ModifiableStructuralVariantCluster segment = null;
        for (StructuralVariantPosition variant : variants) {
            while (ratioIndex < ratios.size() - 1 && ratios.get(ratioIndex).position() < variant.position()) {
                ratioIndex++;
            }

            final long start = start(variant.position(), ratioIndex, ratios);
            final long end = end(variant.position(), ratioIndex, ratios);

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
        assert (index <= ratios.size());
        final long min = position / windowSize * windowSize + 1 - windowSize;
        if (!ratios.isEmpty()) {
            for (int i = index; i >= 0; i--) {
                final CobaltRatio ratio = ratios.get(i);
                if (ratio.position() <= min && Doubles.greaterThan(ratio.tumorGCRatio(), -1)) {
                    return ratio.position();
                }
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

    private static Multimap<String, StructuralVariantPosition> asMap(@NotNull final List<StructuralVariantPosition> variants) {
        final Multimap<String, StructuralVariantPosition> result = ArrayListMultimap.create();

        for (StructuralVariantPosition variant : variants) {
            result.put(variant.chromosome(), variant);
        }

        return result;
    }
}
