package com.hartwig.hmftools.purple.segment;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.window.Window;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.collection.Multimaps;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.common.sv.StructuralVariant;

public class ClusterFactory
{
    private final int mWindowSize;
    private final Window mWindow;

    ClusterFactory(final int windowSize)
    {
        mWindowSize = windowSize;
        mWindow = new Window(windowSize);
    }

    public ListMultimap<Chromosome, Cluster> cluster(final List<StructuralVariant> variants,
            final Multimap<Chromosome, PCFPosition> pcfPositions, final Map<Chromosome,List<CobaltRatio>> ratios)
    {
        final Multimap<Chromosome, SVSegment> positions = Multimaps.fromPositions(SVSegmentFactory.create(variants));
        return cluster(positions, pcfPositions, ratios);
    }

    private ListMultimap<Chromosome, Cluster> cluster(final Multimap<Chromosome, SVSegment> variantPositions,
            final Multimap<Chromosome, PCFPosition> pcfPositions, final Map<Chromosome,List<CobaltRatio>> ratios)
    {
        ListMultimap<Chromosome, Cluster> clusters = ArrayListMultimap.create();
        for(Chromosome chromosome : pcfPositions.keySet())
        {
            final Collection<PCFPosition> chromosomePcfPositions = pcfPositions.get(chromosome);
            final List<CobaltRatio> chromosomeRatios = ratios.containsKey(chromosome) ? ratios.get(chromosome) : Lists.newArrayList();
            final Collection<SVSegment> chromosomeVariants =
                    variantPositions.containsKey(chromosome) ? variantPositions.get(chromosome) : Lists.newArrayList();
            clusters.putAll(chromosome, cluster(chromosomeVariants, chromosomePcfPositions, chromosomeRatios));
        }

        return clusters;
    }

    @VisibleForTesting
    List<Cluster> cluster(
            final Collection<SVSegment> variantPositions, final Collection<PCFPosition> pcfPositions, final List<CobaltRatio> cobaltRatios)
    {
        final List<GenomePosition> allPositions = Lists.newArrayList();
        allPositions.addAll(variantPositions);
        allPositions.addAll(pcfPositions);
        Collections.sort(allPositions);

        final List<Cluster> result = Lists.newArrayList();

        int cobaltIndex = 0;
        Cluster segment = null;
        for(GenomePosition position : allPositions)
        {
            if(position.position() == 1)
            {
                continue;
            }

            while(cobaltIndex < cobaltRatios.size() - 1 && cobaltRatios.get(cobaltIndex).position() < position.position())
            {
                cobaltIndex++;
            }

            final int earliestDetectableCopyNumberChangePosition =
                    earliestDetectableCopyNumberChangePosition(position.position(), cobaltIndex, cobaltRatios);
            if(segment == null || earliestDetectableCopyNumberChangePosition > segment.end())
            {
                if(segment != null)
                {
                    result.add(segment);
                }

                segment = new Cluster(position.chromosome(), earliestDetectableCopyNumberChangePosition, 0);
            }

            segment.End = position.position();

            if(position instanceof SVSegment)
            {
                segment.Variants.add((SVSegment) position);
            }
            else
            {
                segment.PcfPositions.add((PCFPosition) position);
            }
        }
        if(segment != null)
        {
            result.add(segment);
        }

        return result;
    }

    @VisibleForTesting
    int earliestDetectableCopyNumberChangePosition(int position, int index, final List<CobaltRatio> ratios)
    {
        assert (index <= ratios.size());
        final int min = mWindow.start(position) - mWindowSize + 1;
        if(!ratios.isEmpty())
        {
            for(int i = index; i >= 0; i--)
            {
                final CobaltRatio ratio = ratios.get(i);
                if(ratio.position() <= min && Doubles.greaterThan(ratio.tumorGCRatio(), -1))
                {
                    return ratio.position() + 1;
                }
            }
        }

        return min;
    }
}
