package com.hartwig.hmftools.purple.segment;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

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

    public List<PurpleSupportSegment> createSegments(
            final List<StructuralVariant> variants, final Map<Chromosome,List<PCFPosition>> pcfPositions,
            final Map<Chromosome,List<CobaltRatio>> ratios)
    {
        ClusterFactory clusterFactory = new ClusterFactory(mWindowSize);
        final Multimap<Chromosome, Cluster> clusterMap = clusterFactory.cluster(variants, pcfPositions, ratios);
        return segmentCluster(clusterMap);
    }

    private List<PurpleSupportSegment> segmentCluster(final Multimap<Chromosome, Cluster> clusters)
    {
        final List<PurpleSupportSegment> results = Lists.newArrayList();
        results.addAll(segmentMap(clusters).values());
        Collections.sort(results);
        return results;
    }

    private Multimap<Chromosome, PurpleSupportSegment> segmentMap(final Multimap<Chromosome, Cluster> clusters)
    {
        final Multimap<Chromosome, PurpleSupportSegment> segments = ArrayListMultimap.create();
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
    static List<PurpleSupportSegment> create(final GenomePosition centromere, final GenomePosition length, final Collection<Cluster> clusters)
    {
        final List<PurpleSupportSegment> segments = create(length, clusters);

        for(int i = 1; i < segments.size(); ++i)
        {
            PurpleSupportSegment prevSegment = segments.get(i - 1);
            PurpleSupportSegment segment = segments.get(i);

            segment.MinStart = max(segment.MinStart, prevSegment.End + 1);
        }

        return addCentromere(centromere, segments);
    }

    private static List<PurpleSupportSegment> create(final GenomePosition length, final Collection<Cluster> clusters)
    {
        final List<PurpleSupportSegment> result = Lists.newArrayList();
        PurpleSupportSegment segment = create(length.chromosome());
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

    private static PurpleSupportSegment create(String chromosome)
    {
        return new PurpleSupportSegment(chromosome, 1, 0, true, SegmentSupport.NONE, false, 1, 1);
    }

    private static PurpleSupportSegment create(String chromosome, int start, final List<PCFPosition> pcfPositions)
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

        return new PurpleSupportSegment(chromosome, start, 0, true, SegmentSupport.NONE, false, minStart, maxStart);
    }

    private static PurpleSupportSegment createFromCluster(Cluster cluster, SVSegment variant, boolean ratioSupport)
    {
        return new PurpleSupportSegment(
                cluster.chromosome(), variant.Position, 0, ratioSupport, SegmentSupport.fromVariant(variant.Type),
                true, variant.Position, variant.Position);
    }

    private static List<PurpleSupportSegment> addCentromere(final GenomePosition centromere, final List<PurpleSupportSegment> segments)
    {
        final List<PurpleSupportSegment> result = Lists.newArrayList();

        for(PurpleSupportSegment segment : segments)
        {
            if(centromere != null && segment.contains(centromere))
            {
                if(segment.start() == centromere.position())
                {
                    final PurpleSupportSegment start = PurpleSupportSegment.from(segment);
                    start.Support = SegmentSupport.CENTROMERE;
                    result.add(start);
                }
                else
                {
                    final PurpleSupportSegment start = PurpleSupportSegment.from(segment);
                    start.End = centromere.position() - 1;
                    start.MaxStart = min(start.MaxStart, start.End);

                    final PurpleSupportSegment end = PurpleSupportSegment.from(segment);
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

    public static boolean validateSegments(final List<PurpleSupportSegment> segments)
    {
        boolean isValid = true;

        for(int i = 1; i < segments.size(); ++i)
        {
            PurpleSupportSegment segment = segments.get(i);

            if(!positionsWithin(segment.MinStart, segment.MaxStart, segment.start(), segment.end()))
            {
                PPL_LOGGER.error("purple segment({}:{}-{}) has invalid min/maxStart({}-{})",
                        segment.chromosome(), segment.start(), segment.end(),
                        segment.MinStart, segment.MaxStart);

                isValid = false;
            }

            PurpleSupportSegment prevSegment = segments.get(i - 1);

            if(!segment.chromosome().equals(prevSegment.chromosome()))
                continue;

            if(segment.start() <= prevSegment.end())
            {
                PPL_LOGGER.error("purple segment({}:{}-{}) overlaps previous({}-{})",
                        segment.chromosome(), segment.start(), segment.end(), prevSegment.start(), prevSegment.end());
                isValid = false;
            }
        }

        return isValid;
    }
}
