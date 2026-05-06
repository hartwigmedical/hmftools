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
import com.hartwig.hmftools.common.genome.region.Window;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.Multimaps;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.utils.pcf.PCFSource;

public class ClusterFactory
{
    private final int mWindowSize;
    private final Window mWindow;

    ClusterFactory(final int windowSize)
    {
        mWindowSize = windowSize;
        mWindow = new Window(windowSize);
    }

    public ListMultimap<Chromosome, Cluster> cluster(
            final List<StructuralVariant> variants, Map<Chromosome,List<PCFPosition>> pcfPositions,
            final Map<Chromosome,List<CobaltRatio>> ratios)
    {
        final Multimap<Chromosome, SvPosition> positions = Multimaps.fromPositions(SVSegmentFactory.create(variants));
        return buildClusters(positions, pcfPositions, ratios);
    }

    private ListMultimap<Chromosome, Cluster> buildClusters(
            final Multimap<Chromosome, SvPosition> variantPositions, final Map<Chromosome,List<PCFPosition>> chrPcfPositions,
            final Map<Chromosome,List<CobaltRatio>> ratios)
    {
        ListMultimap<Chromosome, Cluster> chrClusters = ArrayListMultimap.create();

        for(Chromosome chromosome : chrPcfPositions.keySet())
        {
            List<PCFPosition> pcfPositions = chrPcfPositions.get(chromosome);
            List<CobaltRatio> cobaltRatios = ratios.containsKey(chromosome) ? ratios.get(chromosome) : Collections.emptyList();

            Collection<SvPosition> svPositions = variantPositions.containsKey(chromosome)
                    ? variantPositions.get(chromosome) : Collections.emptyList();

            List<Cluster> clusters = buildChromosomeClusters(svPositions, pcfPositions, cobaltRatios);
            chrClusters.putAll(chromosome, clusters);
        }

        return chrClusters;
    }

    @VisibleForTesting
    List<Cluster> buildChromosomeClusters(
            final Collection<SvPosition> svPositions, final List<PCFPosition> pcfPositions, final List<CobaltRatio> cobaltRatios)
    {
        List<GenomePosition> allPositions = Lists.newArrayList();
        allPositions.addAll(svPositions);
        allPositions.addAll(pcfPositions);
        Collections.sort(allPositions);

        List<Cluster> clusters = Lists.newArrayList();

        int cobaltIndex = 0;
        Cluster segment = null;
        PCFPosition lastAmberPcf = null;
        GenomePosition lastPosition = null;

        for(GenomePosition position : allPositions)
        {
            if(position.position() == 1)
                continue;

            while(cobaltIndex < cobaltRatios.size() - 1 && cobaltRatios.get(cobaltIndex).position() < position.position())
            {
                cobaltIndex++;
            }

            boolean canSegment = true;

            if(position instanceof PCFPosition)
            {
                PCFPosition pcfPosition = (PCFPosition)position;

                if(pcfPosition.Source == PCFSource.TUMOR_BAF)
                {
                    if(lastAmberPcf != null && lastPosition != null && lastPosition.position() > lastAmberPcf.position())
                        canSegment = false;

                    lastAmberPcf = pcfPosition;
                }
            }

            int earliestCnChangePosition = 0;
            if(canSegment)
            {
                earliestCnChangePosition = findFirstValidCobaltRatio(position.position(), cobaltIndex, cobaltRatios);

                if(segment != null)
                    canSegment = earliestCnChangePosition > segment.end();
            }

            if(segment == null || canSegment)
            {
                segment = new Cluster(position.chromosome(), earliestCnChangePosition, 0);
                clusters.add(segment);
            }

            segment.End = position.position();

            if(position instanceof SvPosition)
            {
                segment.Variants.add((SvPosition) position);
            }
            else
            {
                segment.PcfPositions.add((PCFPosition)position);
            }

            lastPosition = position;
        }

        return clusters;
    }

    @VisibleForTesting
    int findFirstValidCobaltRatio(int position, int index, final List<CobaltRatio> ratios)
    {
        // returns the first valid ratio earlier than the specific position, starting at the specified index and working backwards
        int min = mWindow.start(position) - mWindowSize + 1;

        if(!ratios.isEmpty())
        {
            for(int i = index; i >= 0; i--)
            {
                CobaltRatio ratio = ratios.get(i);
                if(ratio.position() <= min && Doubles.greaterThan(ratio.tumorGCRatio(), -1))
                {
                    return ratio.position() + 1;
                }
            }
        }

        return min;
    }
}
