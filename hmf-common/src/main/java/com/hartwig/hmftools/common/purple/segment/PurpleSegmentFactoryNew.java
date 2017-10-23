package com.hartwig.hmftools.common.purple.segment;

import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.ChromosomeLength;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

public class PurpleSegmentFactoryNew {

    @NotNull
    public static List<PurpleSegment> segment(@NotNull final Multimap<String, StructuralVariantCluster> clusters,
            @NotNull final Multimap<String, GenomePosition> ratios, @NotNull final Map<String, ChromosomeLength> lengths) {

        final List<PurpleSegment> segments = Lists.newArrayList();

        for (String chromosome : lengths.keySet()) {
            if (HumanChromosome.contains(chromosome)) {
                final Collection<StructuralVariantCluster> cluster =
                        clusters.containsKey(chromosome) ? clusters.get(chromosome) : Collections.emptyList();
                final Collection<GenomePosition> ratio = ratios.containsKey(chromosome) ? ratios.get(chromosome) : Collections.emptyList();
                segments.addAll(create(lengths.get(chromosome), cluster, ratio));
            }
        }

        Collections.sort(segments);
        return segments;
    }

    @NotNull
    public static List<PurpleSegment> create(@NotNull final ChromosomeLength chromosome,
            @NotNull final Collection<StructuralVariantCluster> clusteredVariants, @NotNull final Collection<GenomePosition> ratioBreaks) {

        final List<PurpleSegment> result = Lists.newArrayList();

        Iterator<GenomePosition> ratioIterator = ratioBreaks.iterator();
        Iterator<StructuralVariantCluster> clusterIterator = clusteredVariants.iterator();

        GenomePosition ratio = ratioIterator.hasNext() ? ratioIterator.next() : null;
        StructuralVariantCluster cluster = clusterIterator.hasNext() ? clusterIterator.next() : null;

        ModifiablePurpleSegment segment = create(chromosome.chromosome(), 1);
        while (ratio != null || cluster != null) {

            if (ratio == null) {
                result.add(segment.setEnd(cluster.firstVariantPosition() - 1));
                segment = create(cluster);

                if (cluster.firstVariantPosition() != cluster.finalVariantPosition()) {
                    result.add(segment);
                    segment = createFromLastVariant(cluster, false);
                }
                cluster = clusterIterator.hasNext() ? clusterIterator.next() : null;

            } else if (cluster == null) {
                result.add(segment.setEnd(ratio.position() - 1));
                segment = create(ratio.chromosome(), ratio.position());
                ratio = ratioIterator.hasNext() ? ratioIterator.next() : null;
            } else {

                if (cluster.start() <= ratio.position()) {
                    result.add(segment.setEnd(cluster.firstVariantPosition() - 1));
                    segment = create(cluster);
                    final boolean ratioSupport;
                    if (ratio.position() <= segment.maxBoundary()) {
                        ratioSupport = true;
                        segment.setRatioSupport(true);
                        ratio = ratioIterator.hasNext() ? ratioIterator.next() : null;
                    } else {
                        ratioSupport = false;
                    }

                    if (cluster.firstVariantPosition() != cluster.finalVariantPosition()) {
                        result.add(segment);
                        segment = createFromLastVariant(cluster, ratioSupport);
                    }
                    cluster = clusterIterator.hasNext() ? clusterIterator.next() : null;
                } else {
                    if (ratio.position() <= segment.maxBoundary()) {
                        segment.setRatioSupport(true);
                    } else {
                        result.add(segment.setEnd(ratio.position() - 1));
                        segment = create(ratio.chromosome(), ratio.position());
                    }
                    ratio = ratioIterator.hasNext() ? ratioIterator.next() : null;
                }
            }
        }

        result.add(segment.setEnd(chromosome.position()));
        return result;
    }

    private static ModifiablePurpleSegment create(String chromosome, long start) {
        return ModifiablePurpleSegment.create()
                .setChromosome(chromosome)
                .setRatioSupport(true)
                .setStart(start)
                .setEnd(0)
                .setMaxBoundary(0)
                .setStructuralVariantSupport(StructuralVariantSupport.NONE);
    }

    private static ModifiablePurpleSegment create(StructuralVariantCluster cluster) {
        return ModifiablePurpleSegment.create()
                .setChromosome(cluster.chromosome())
                .setRatioSupport(false)
                .setStart(cluster.firstVariantPosition())
                .setEnd(cluster.finalVariantPosition() - 1)
                .setMaxBoundary(cluster.end())
                .setStructuralVariantSupport(cluster.type());
    }

    private static ModifiablePurpleSegment createFromLastVariant(StructuralVariantCluster cluster, boolean ratioSupport) {
        return ModifiablePurpleSegment.create()
                .setChromosome(cluster.chromosome())
                .setRatioSupport(ratioSupport)
                .setStart(cluster.finalVariantPosition())
                .setEnd(0)
                .setMaxBoundary(cluster.end())
                .setStructuralVariantSupport(cluster.finalVariantType());
    }
}
