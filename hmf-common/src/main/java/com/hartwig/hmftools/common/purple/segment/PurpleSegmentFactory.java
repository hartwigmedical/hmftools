package com.hartwig.hmftools.common.purple.segment;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.pcf.PCFPosition;
import com.hartwig.hmftools.common.pcf.PCFSource;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class PurpleSegmentFactory {

    private final int windowSize;
    private final Map<Chromosome, GenomePosition> lengths;
    private final Map<Chromosome, GenomePosition> centromeres;

    public PurpleSegmentFactory(final int windowSize, final Map<Chromosome, GenomePosition> centromeres, final Map<Chromosome, GenomePosition> lengths) {
        this.windowSize = windowSize;
        this.centromeres = centromeres;
        this.lengths = lengths;
    }

    public List<PurpleSegment> segment(@NotNull final List<StructuralVariant> variants,
            @NotNull final Multimap<Chromosome, PCFPosition> pcfPositions, @NotNull ListMultimap<Chromosome, CobaltRatio> ratios) {
        final Multimap<Chromosome, Cluster> clusterMap = new ClusterFactory(windowSize).cluster(variants, pcfPositions, ratios);
        return segmentCluster(clusterMap);
    }

    @NotNull
    private List<PurpleSegment> segmentCluster(@NotNull final Multimap<Chromosome, Cluster> clusters) {
        final List<PurpleSegment> results = Lists.newArrayList();
        results.addAll(segmentMap(clusters).values());
        Collections.sort(results);
        return results;
    }

    @NotNull
    private Multimap<Chromosome, PurpleSegment> segmentMap(@NotNull final Multimap<Chromosome, Cluster> clusters) {

        final Multimap<Chromosome, PurpleSegment> segments = ArrayListMultimap.create();
        for (Chromosome chromosome : clusters.keySet()) {

            GenomePosition length = lengths.get(chromosome);
            GenomePosition centromere = centromeres.get(chromosome);

            final Collection<Cluster> cluster = clusters.containsKey(chromosome) ? clusters.get(chromosome) : Collections.emptyList();
            segments.putAll(chromosome, create(centromere, length, cluster));
        }

        return segments;
    }

    @NotNull
    @VisibleForTesting
    static List<PurpleSegment> create(@NotNull final GenomePosition centromere, @NotNull final GenomePosition length,
            @NotNull final Collection<Cluster> clusters) {
        final List<PurpleSegment> result = Lists.newArrayList();
        ModifiablePurpleSegment segment = create(length.chromosome(), 1).setSupport(SegmentSupport.TELOMERE);

        for (final Cluster cluster : clusters) {
            boolean ratioSupport = !cluster.ratios().isEmpty();

            final List<ClusterVariantLeg> variants = cluster.variants();
            if (!variants.isEmpty()) {
                for (final ClusterVariantLeg variant : variants) {
                    if (variant.position() != segment.start()) {
                        result.add(setEnd(centromere, segment, (variant.position() - 1)));
                        segment = createFromCluster(cluster, variant, ratioSupport);
                    } else {
                        segment.setSupport(SegmentSupport.MULTIPLE);
                    }
                }
                segment.setSvCluster(false);
            } else {

                final List<PCFPosition> pcfPositions = cluster.pcfPositions();

                // JOBA: DO FIRST
                final GenomePosition firstRatioBreak = pcfPositions.get(0);
                result.add(setEnd(centromere, segment, firstRatioBreak.position() - 1));
                segment = create(firstRatioBreak.chromosome(), firstRatioBreak.position(), pcfPositions);
            }
        }

        result.add(setEnd(centromere, segment, length.position()));
        return result;
    }

    @NotNull
    private static ModifiablePurpleSegment create(@NotNull String chromosome, long start) {
        return ModifiablePurpleSegment.create()
                .setChromosome(chromosome)
                .setRatioSupport(true)
                .setStart(start)
                .setMinStart(start)
                .setMaxStart(start)
                .setEnd(0)
                .setSvCluster(false)
                .setSupport(SegmentSupport.NONE);
    }

    @NotNull
    private static ModifiablePurpleSegment create(@NotNull String chromosome, long start, @NotNull final List<PCFPosition> pcfPositions) {

        long minStart = pcfPositions.stream()
                .filter(x -> x.source() == PCFSource.TUMOR_RATIO)
                .mapToLong(PCFPosition::minPosition)
                .min()
                .orElse(start);
        long maxStart = pcfPositions.stream()
                .filter(x -> x.source() == PCFSource.TUMOR_RATIO)
                .mapToLong(PCFPosition::maxPosition)
                .max()
                .orElse(start);

        return ModifiablePurpleSegment.create()
                .setChromosome(chromosome)
                .setRatioSupport(true)
                .setStart(start)
                .setMinStart(minStart)
                .setMaxStart(maxStart)
                .setEnd(0)
                .setSvCluster(false)
                .setSupport(SegmentSupport.NONE);
    }

    @NotNull
    private static ModifiablePurpleSegment createFromCluster(@NotNull Cluster cluster, @NotNull ClusterVariantLeg variant,
            boolean ratioSupport) {
        return ModifiablePurpleSegment.create()
                .setChromosome(cluster.chromosome())
                .setRatioSupport(ratioSupport)
                .setStart(variant.position())
                .setMinStart(variant.position())
                .setMaxStart(variant.position())
                .setEnd(0)
                .setSvCluster(true)
                .setSupport(SegmentSupport.fromVariant(variant.type()));
    }

    private static ModifiablePurpleSegment setEnd(@Nullable GenomePosition centromere, @NotNull ModifiablePurpleSegment segment, long end) {
        segment.setEnd(end);

        if (centromere != null && segment.contains(centromere)) {
            segment.setSupport(SegmentSupport.CENTROMERE);
        }

        return segment;
    }
}
