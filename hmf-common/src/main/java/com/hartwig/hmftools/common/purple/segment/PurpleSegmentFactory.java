package com.hartwig.hmftools.common.purple.segment;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.centromeres.Centromeres;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.ChromosomeLength;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.pcf.PCFPosition;
import com.hartwig.hmftools.common.pcf.PCFSource;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public final class PurpleSegmentFactory {

    private static final Map<String, GenomeRegion> CENTROMERES = Centromeres.grch37();

    @NotNull
    public static List<PurpleSegment> segment(@NotNull final Multimap<String, Cluster> clusters,
            @NotNull final Map<Chromosome, ChromosomeLength> lengths) {
        final List<PurpleSegment> results = Lists.newArrayList();
        results.addAll(segmentMap(clusters, lengths).values());
        Collections.sort(results);
        return results;
    }

    @NotNull
    private static Multimap<String, PurpleSegment> segmentMap(@NotNull final Multimap<String, Cluster> clusters,
            @NotNull final Map<Chromosome, ChromosomeLength> lengths) {
        final Multimap<String, PurpleSegment> segments = ArrayListMultimap.create();
        for (String chromosome : clusters.keySet()) {
            if (HumanChromosome.contains(chromosome)) {
                final Collection<Cluster> cluster = clusters.containsKey(chromosome) ? clusters.get(chromosome) : Collections.emptyList();
                segments.putAll(chromosome, create(lengths.get(HumanChromosome.fromString(chromosome)), cluster));
            }
        }

        return segments;
    }

    @NotNull
    @VisibleForTesting
    static List<PurpleSegment> create(@NotNull final ChromosomeLength chromosome, @NotNull final Collection<Cluster> clusters) {
        final List<PurpleSegment> result = Lists.newArrayList();
        ModifiablePurpleSegment segment = create(chromosome.chromosome(), 1).setSupport(SegmentSupport.TELOMERE);

        for (final Cluster cluster : clusters) {
            boolean ratioSupport = !cluster.ratios().isEmpty();

            final List<ClusterVariantLeg> variants = cluster.variants();
            if (!variants.isEmpty()) {
                for (final ClusterVariantLeg variant : variants) {
                    if (variant.position() != segment.start()) {
                        result.add(setEnd(segment, (variant.position() - 1)));
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
                result.add(setEnd(segment, firstRatioBreak.position() - 1));
                segment = create(firstRatioBreak.chromosome(), firstRatioBreak.position(), pcfPositions);
            }
        }

        result.add(setEnd(segment, chromosome.length()));
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

    private static ModifiablePurpleSegment setEnd(@NotNull ModifiablePurpleSegment segment, long end) {
        segment.setEnd(end);

        final GenomeRegion centromere = CENTROMERES.get(segment.chromosome());
        if (centromere != null && centromere.overlaps(segment)) {
            segment.setSupport(SegmentSupport.CENTROMERE);
        }

        return segment;
    }
}
