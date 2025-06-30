package com.hartwig.hmftools.linx.visualiser.circos;

import static com.hartwig.hmftools.linx.visualiser.SvVisualiser.VIS_LOGGER;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

public class Downsampler
{
    private static final int DEFAULT_DOWNSAMPLE_TARGET = 5000;
    private static final int DEFAULT_MIN_PRESERVED = 5;

    public static <T> List<T> downsampleList(List<T> list, int downsampleTarget)
    {
        if(list.size() <= downsampleTarget)
            return list;

        double step = (double) list.size() / downsampleTarget;

        List<T> subsampledList = Lists.newArrayList();

        for (int sampleIndex = 0; sampleIndex < downsampleTarget; sampleIndex++) {
            int newSampleIndex = (int) Math.floor(sampleIndex * step);

            if (newSampleIndex >= list.size())
                break;

            subsampledList.add(list.get(newSampleIndex));
        }

        return subsampledList;
    }

    public static int[] calcGroupDownsampleTargets(List<Integer> groupSizes, int downsampleTarget, int minPreserved)
    {
        int nGroups = groupSizes.size();

        int[] reservedCounts = new int[nGroups];
        int[] extraCounts = new int[nGroups];

        for(int i = 0; i < nGroups; i++)
        {
            int groupSize = groupSizes.get(i);

            if(groupSize > 0)
            {
                reservedCounts[i] = minPreserved;
                extraCounts[i] = groupSize - minPreserved;
            }
        }

        int reservedSum = Arrays.stream(reservedCounts).sum();
        int extraSum = Arrays.stream(extraCounts).sum();

        int countsToDistribute = downsampleTarget - reservedSum;

        int[] extraCountsDownsampled = new int[nGroups];
        for(int i = 0; i < nGroups; i++)
        {
            int extraCount = extraCounts[i];
            double extraFraction = (double) extraCount / extraSum;

            extraCountsDownsampled[i] = (int) Math.round(extraFraction * countsToDistribute);
            // This^ can cause rounding errors, but we don't need to be precise with the amount to downsample
        }

        int[] groupSizesDownsampled = new int[nGroups];
        for(int i = 0; i < nGroups; i++)
        {
            groupSizesDownsampled[i] = reservedCounts[i] + extraCountsDownsampled[i];
        }

        return groupSizesDownsampled;
    }

    public static List<GenomeRegion> regionsFromPositions(List<GenomePosition> positionsToScale)
    {
        List<String> chromosomes = positionsToScale.stream().map(x -> x.chromosome()).distinct().toList();

        List<GenomeRegion> contigs = Lists.newArrayList();

        for(String chromosome : chromosomes)
        {
            List<GenomePosition> positions = positionsToScale.stream()
                    .filter(x -> x.chromosome().equals(chromosome))
                    .collect(Collectors.toList());

            positions.add(GenomePositions.create(chromosome, 1));
            positions.add(GenomePositions.create(chromosome, Integer.MAX_VALUE));

            positions = positions.stream().distinct().sorted().toList();

            for(int i = 1; i < positions.size(); i++)
            {
                int start = positions.get(i-1).position();
                if(start == 1) start = 0; // Convert to 0-based coordinates

                int end = positions.get(i).position();
                GenomeRegion contig = GenomeRegions.create(chromosome, start, end);
                contigs.add(contig);
            }
        }

        return contigs;
    }

    public static <T extends GenomePosition> List<T> downsampleWithMinimumPerContig(List<T> dataPoints, List<GenomePosition> positionsToScale)
    {
        if(dataPoints.size() < DEFAULT_DOWNSAMPLE_TARGET)
            return dataPoints;

        List<GenomeRegion> contigs = regionsFromPositions(positionsToScale);

        List<List<T>> dataPointsGrouped = Lists.newArrayList();
        for(GenomeRegion contig : contigs)
        {
            List<T> contigAmberBAFs = dataPoints.stream()
                    .filter(x ->
                            x.chromosome().equals(contig.chromosome()) &&
                            x.position() > contig.start() && x.position() <= contig.end()
                    ).toList();

            dataPointsGrouped.add(contigAmberBAFs);
        }

        List<Integer> groupSizes = dataPointsGrouped.stream().map(x -> x.size()).toList();
        int[] groupDownsampleTargets = calcGroupDownsampleTargets(groupSizes, DEFAULT_DOWNSAMPLE_TARGET, DEFAULT_MIN_PRESERVED);

        List<T> result = Lists.newArrayList();
        for(int i = 0; i < contigs.size(); i++)
        {
            List<T> groupDataPoints = dataPointsGrouped.get(i);
            int groupDownsampleTarget = groupDownsampleTargets[i];

            List<T> groupDataPointsDownsampled = downsampleList(groupDataPoints, groupDownsampleTarget);
            result.addAll(groupDataPointsDownsampled);

            boolean isDownsampled = groupDataPoints.size() > groupDownsampleTarget;
            if(isDownsampled)
            {
                VIS_LOGGER.trace("region {}:{}-{} downsampled {} -> {}, first data point: {}",
                        contigs.get(i).chromosome(), contigs.get(i).start(), contigs.get(i).end(),
                        groupDataPoints.size(), groupDataPointsDownsampled.size(),
                        dataPoints.get(0)
                );
            }
        }

        return result;
    }

}
