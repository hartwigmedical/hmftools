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
import com.hartwig.hmftools.common.purple.SegmentSupport;
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
        final List<PurpleSegment> segments = create(length, clusters);
        return addCentromere(centromere, segments);
    }

    private static List<PurpleSegment> create(final GenomePosition length, final Collection<Cluster> clusters)
    {
        final List<PurpleSegment> result = Lists.newArrayList();
        PurpleSegment segment = create(length.chromosome());
        segment.Support = SegmentSupport.TELOMERE;

        for(final Cluster cluster : clusters)
        {
            boolean ratioSupport = !cluster.ratios().isEmpty();

            final List<SVSegment> variants = cluster.Variants;
            if(!variants.isEmpty())
            {
                for(final SVSegment variant : variants)
                {
                    if(variant.position() != segment.start())
                    {
                        segment.End = variant.position() - 1;
                        result.add(segment);
                        segment = createFromCluster(cluster, variant, ratioSupport);
                    }
                    else
                    {
                        segment.Support = SegmentSupport.MULTIPLE;
                    }
                }

                segment.SvCluster = false;
            }
            else
            {

                final List<PCFPosition> pcfPositions = cluster.PcfPositions;

                // DO FIRST
                final GenomePosition firstRatioBreak = pcfPositions.get(0);
                segment.End = firstRatioBreak.position() - 1;
                result.add(segment);
                segment = create(firstRatioBreak.chromosome(), firstRatioBreak.position(), pcfPositions);
            }
        }

        segment.End = length.position();
        result.add(segment);
        return result;
    }

    private static PurpleSegment create(String chromosome)
    {
        return new PurpleSegment(chromosome, 1, 0, true, SegmentSupport.NONE, false, 1, 1);
    }

    private static PurpleSegment create(String chromosome, int start, final List<PCFPosition> pcfPositions)
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

        return new PurpleSegment(chromosome, start, 0, true, SegmentSupport.NONE, false, minStart, maxStart);
    }

    private static PurpleSegment createFromCluster(Cluster cluster, SVSegment variant, boolean ratioSupport)
    {
        return new PurpleSegment(
                cluster.chromosome(), variant.Position, 0, ratioSupport, SegmentSupport.fromVariant(variant.Type),
                true, variant.Position, variant.Position);
    }

    private static List<PurpleSegment> addCentromere(final GenomePosition centromere, final List<PurpleSegment> segments)
    {
        final List<PurpleSegment> result = Lists.newArrayList();

        for(PurpleSegment segment : segments)
        {
            if(centromere != null && segment.contains(centromere))
            {
                if(segment.start() == centromere.position())
                {
                    final PurpleSegment start = PurpleSegment.from(segment);
                    start.Support = SegmentSupport.CENTROMERE;
                    result.add(start);
                }
                else
                {
                    final PurpleSegment start = PurpleSegment.from(segment);
                    start.End = centromere.position() - 1;

                    final PurpleSegment end = PurpleSegment.from(segment);
                    end.Support = SegmentSupport.CENTROMERE;
                    end.Start = centromere.position();
                    end.MinStart = centromere.position();
                    end.MaxStart = centromere.position();

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
