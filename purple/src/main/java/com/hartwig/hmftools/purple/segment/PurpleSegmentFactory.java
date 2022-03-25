package com.hartwig.hmftools.purple.segment;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.common.utils.pcf.PCFSource;
import com.hartwig.hmftools.common.sv.StructuralVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PurpleSegmentFactory
{
    private final int mWindowSize;
    private final Map<Chromosome, GenomePosition> mLengths;
    private final Map<Chromosome, GenomePosition> mCentromeres;

    public PurpleSegmentFactory(final int windowSize, final Map<Chromosome, GenomePosition> centromeres,
            final Map<Chromosome, GenomePosition> lengths)
    {
        mWindowSize = windowSize;
        mCentromeres = centromeres;
        mLengths = lengths;
    }

    public List<PurpleSegment> segment(
            final List<StructuralVariant> variants, final Multimap<Chromosome, PCFPosition> pcfPositions,
            final Map<Chromosome,List<CobaltRatio>> ratios)
    {
        ClusterFactory clusterFactory = new ClusterFactory(mWindowSize);
        final Multimap<Chromosome, Cluster> clusterMap = clusterFactory.cluster(variants, pcfPositions, ratios);
        return segmentCluster(clusterMap);
    }

    @NotNull
    private List<PurpleSegment> segmentCluster(final Multimap<Chromosome, Cluster> clusters)
    {
        final List<PurpleSegment> results = Lists.newArrayList();
        results.addAll(segmentMap(clusters).values());
        Collections.sort(results);
        return results;
    }

    private Multimap<Chromosome, PurpleSegment> segmentMap(final Multimap<Chromosome, Cluster> clusters)
    {
        final Multimap<Chromosome, PurpleSegment> segments = ArrayListMultimap.create();
        for(Chromosome chromosome : clusters.keySet())
        {
            GenomePosition length = mLengths.get(chromosome);
            GenomePosition centromere = mCentromeres.get(chromosome);

            final Collection<Cluster> cluster = clusters.containsKey(chromosome) ? clusters.get(chromosome) : Collections.emptyList();
            segments.putAll(chromosome, create(centromere, length, cluster));
        }

        return segments;
    }

    @VisibleForTesting
    static List<PurpleSegment> create(final GenomePosition centromere, final GenomePosition length, final Collection<Cluster> clusters)
    {
        return addCentromere(centromere, create(length, clusters));
    }

    private static List<PurpleSegment> create(final GenomePosition length, final Collection<Cluster> clusters)
    {
        final List<PurpleSegment> result = Lists.newArrayList();
        ModifiablePurpleSegment segment = create(length.chromosome()).setSupport(SegmentSupport.TELOMERE);

        for(final Cluster cluster : clusters)
        {
            boolean ratioSupport = !cluster.ratios().isEmpty();

            final List<SVSegment> variants = cluster.variants();
            if(!variants.isEmpty())
            {
                for(final SVSegment variant : variants)
                {
                    if(variant.position() != segment.start())
                    {
                        result.add(segment.setEnd(variant.position() - 1));
                        segment = createFromCluster(cluster, variant, ratioSupport);
                    }
                    else
                    {
                        segment.setSupport(SegmentSupport.MULTIPLE);
                    }
                }
                segment.setSvCluster(false);
            }
            else
            {

                final List<PCFPosition> pcfPositions = cluster.pcfPositions();

                // DO FIRST
                final GenomePosition firstRatioBreak = pcfPositions.get(0);
                result.add(segment.setEnd(firstRatioBreak.position() - 1));
                segment = create(firstRatioBreak.chromosome(), firstRatioBreak.position(), pcfPositions);
            }
        }

        result.add(segment.setEnd(length.position()));
        return result;
    }

    @NotNull
    private static ModifiablePurpleSegment create(String chromosome)
    {
        return ModifiablePurpleSegment.create()
                .setChromosome(chromosome)
                .setRatioSupport(true)
                .setStart(1)
                .setMinStart(1)
                .setMaxStart(1)
                .setEnd(0)
                .setSvCluster(false)
                .setSupport(SegmentSupport.NONE);
    }

    @NotNull
    private static ModifiablePurpleSegment create(String chromosome, int start, final List<PCFPosition> pcfPositions)
    {
        int minStart = pcfPositions.stream()
                .filter(x -> x.Source == PCFSource.TUMOR_RATIO)
                .mapToInt(PCFPosition::minPosition)
                .min()
                .orElse(start);
        int maxStart = pcfPositions.stream()
                .filter(x -> x.Source == PCFSource.TUMOR_RATIO)
                .mapToInt(PCFPosition::maxPosition)
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
    private static ModifiablePurpleSegment createFromCluster(Cluster cluster, SVSegment variant, boolean ratioSupport)
    {
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

    @NotNull
    private static List<PurpleSegment> addCentromere(@Nullable final GenomePosition centromere,
            final List<PurpleSegment> segments)
    {
        final List<PurpleSegment> result = Lists.newArrayList();

        for(PurpleSegment segment : segments)
        {
            if(centromere != null && segment.contains(centromere))
            {
                if(segment.start() == centromere.position())
                {
                    final PurpleSegment start = ImmutablePurpleSegment.builder().from(segment).support(SegmentSupport.CENTROMERE).build();
                    result.add(start);
                }
                else
                {
                    final PurpleSegment start = ImmutablePurpleSegment.builder().from(segment).end(centromere.position() - 1).build();
                    final PurpleSegment end = ImmutablePurpleSegment.builder()
                            .from(segment)
                            .start(centromere.position())
                            .minStart(centromere.position())
                            .maxStart(centromere.position())
                            .support(SegmentSupport.CENTROMERE)
                            .build();

                    result.add(start);
                    result.add(end);
                }

            }
            else
            {
                result.add(segment);
            }

        }

        return result;
    }
}
