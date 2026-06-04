package com.hartwig.hmftools.purple.segment;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

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

public class PurpleSupportSegmentFactory
{
    private final int mWindowSize;
    private final Map<Chromosome, GenomePosition> mLengths;
    private final Map<Chromosome, GenomePosition> mCentromeres;

    public PurpleSupportSegmentFactory(
            final int windowSize, final Map<Chromosome, GenomePosition> centromeres, final Map<Chromosome,GenomePosition> lengths)
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
        Multimap<Chromosome, Cluster> clusterMap = clusterFactory.cluster(variants, pcfPositions, ratios);

        Multimap<Chromosome, PurpleSupportSegment> segments = ArrayListMultimap.create();
        for(Chromosome chromosome : clusterMap.keySet())
        {
            GenomePosition length = mLengths.get(chromosome);
            GenomePosition centromere = mCentromeres.get(chromosome);

            Collection<Cluster> cluster = clusterMap.containsKey(chromosome) ? clusterMap.get(chromosome) : Collections.emptyList();

            List<PurpleSupportSegment> chrSegments = createChromosomeSegments(centromere, length, cluster);
            segments.putAll(chromosome, chrSegments);
        }

        List<PurpleSupportSegment> supportSegments = Lists.newArrayList();
        segments.values().stream().forEach(x -> supportSegments.add(x));
        Collections.sort(supportSegments);

        return supportSegments;
    }

    @VisibleForTesting
    static List<PurpleSupportSegment> createChromosomeSegments(
            final GenomePosition centromere, final GenomePosition length, final Collection<Cluster> clusters)
    {
        List<PurpleSupportSegment> segments = createClusterSegments(length, clusters);

        for(int i = 1; i < segments.size(); ++i)
        {
            PurpleSupportSegment prevSegment = segments.get(i - 1);
            PurpleSupportSegment segment = segments.get(i);

            segment.setMinStart(max(segment.minStart(), prevSegment.end() + 1));
        }

        return addCentromere(centromere, segments);
    }

    private static List<PurpleSupportSegment> createClusterSegments(final GenomePosition length, final Collection<Cluster> clusters)
    {
        List<PurpleSupportSegment> result = Lists.newArrayList();

        PurpleSupportSegment segment = new PurpleSupportSegment(
                length.chromosome(), 1, 0, true, SegmentSupport.NONE, false, 1, 1);

        segment.Support = SegmentSupport.TELOMERE;

        for(Cluster cluster : clusters)
        {
            boolean ratioSupport = !cluster.ratios().isEmpty();

            List<SvPosition> variants = cluster.Variants;
            if(!variants.isEmpty())
            {
                for(SvPosition variant : variants)
                {
                    if(variant.position() != segment.start())
                    {
                        segment.setEnd(variant.position() - 1);
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

                List<PCFPosition> pcfPositions = cluster.PcfPositions;

                // DO FIRST
                GenomePosition firstRatioBreak = pcfPositions.get(0);

                segment.setEnd(firstRatioBreak.position() - 1);

                result.add(segment);
                segment = createPcfSegment(firstRatioBreak.chromosome(), firstRatioBreak.position(), pcfPositions);
            }
        }

        segment.setEnd(length.position());
        result.add(segment);
        return result;
    }

    private static PurpleSupportSegment createPcfSegment(final String chromosome, int start, final List<PCFPosition> pcfPositions)
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

    private static PurpleSupportSegment createFromCluster(Cluster cluster, SvPosition variant, boolean ratioSupport)
    {
        return new PurpleSupportSegment(
                cluster.chromosome(), variant.Position, 0, ratioSupport, SegmentSupport.fromVariant(variant.Type),
                true, variant.Position, variant.Position);
    }

    private static List<PurpleSupportSegment> addCentromere(final GenomePosition centromere, final List<PurpleSupportSegment> segments)
    {
        List<PurpleSupportSegment> result = Lists.newArrayList();

        for(PurpleSupportSegment segment : segments)
        {
            if(centromere != null && segment.contains(centromere))
            {
                if(segment.start() == centromere.position())
                {
                    PurpleSupportSegment start = PurpleSupportSegment.from(segment);
                    start.Support = SegmentSupport.CENTROMERE;
                    result.add(start);
                }
                else
                {
                    PurpleSupportSegment start = PurpleSupportSegment.from(segment);
                    start.setEnd(centromere.position() - 1);
                    start.setMaxStart(min(start.maxStart(), start.end()));

                    PurpleSupportSegment end = PurpleSupportSegment.from(segment);
                    end.Support = SegmentSupport.CENTROMERE;
                    end.setStart(centromere.position());
                    end.setMinStart(centromere.position());
                    end.setMaxStart(centromere.position());

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

            if(!positionsWithin(segment.minStart(), segment.maxStart(), segment.start(), segment.end()))
            {
                PPL_LOGGER.error("purple segment({}:{}-{}) has invalid min/maxStart({}-{})",
                        segment.chromosome(), segment.start(), segment.end(), segment.minStart(), segment.maxStart());

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
