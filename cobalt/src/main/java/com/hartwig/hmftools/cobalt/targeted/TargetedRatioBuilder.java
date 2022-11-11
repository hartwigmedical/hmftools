package com.hartwig.hmftools.cobalt.targeted;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.cobalt.Chromosome;
import com.hartwig.hmftools.cobalt.CobaltConstants;
import com.hartwig.hmftools.cobalt.ratio.RatioBuilder;
import com.hartwig.hmftools.common.cobalt.ImmutableReadRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.window.Window;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.logging.log4j.Level;
import org.jetbrains.annotations.Nullable;

public class TargetedRatioBuilder implements RatioBuilder
{
    private final ArrayListMultimap<Chromosome,ReadRatio> mOnTargetRatios;
    private final ArrayListMultimap<Chromosome,ReadRatio> mOffTargetRatios;
    private final ArrayListMultimap<Chromosome,ReadRatio> mCombinedRatios;

    public TargetedRatioBuilder(
            final List<GenomePosition> targetRegions,
            final Map<GenomePosition, Double> targetEnrichmentRatios, final ListMultimap<Chromosome, ReadRatio> rawRatios)
    {
        mOnTargetRatios = ArrayListMultimap.create();
        mOffTargetRatios = ArrayListMultimap.create();
        mCombinedRatios = ArrayListMultimap.create();

        // calculate median for sample
        populateOnTargetRatios(targetEnrichmentRatios, rawRatios);
        populateOffTargetRatios(rawRatios, CobaltConstants.OFF_TARGET_WINDOW_SIZE, targetRegions);
        CB_LOGGER.info("{} on target GC ratios, {} off target GC ratios", mOnTargetRatios.size(), mOffTargetRatios.size());

        populateCombinedRatios(mOnTargetRatios, mOffTargetRatios);
    }

    // we use on target ratios only for now
    @Override
    public ArrayListMultimap<Chromosome, ReadRatio> ratios() { return mOnTargetRatios; }

    public ArrayListMultimap<Chromosome, ReadRatio> onTargetRatios() { return mOnTargetRatios; }
    public ArrayListMultimap<Chromosome, ReadRatio> offTargetRatios() { return mOffTargetRatios; }

    private void populateCombinedRatios(
            final ArrayListMultimap<Chromosome, ReadRatio> ratios1, final ArrayListMultimap<Chromosome, ReadRatio> ratios2)
    {
        mCombinedRatios.clear();

        Set<Chromosome> chromosomes = Sets.newIdentityHashSet();

        chromosomes.addAll(ratios1.keySet());
        chromosomes.addAll(ratios2.keySet());

        for(Chromosome chromosome : chromosomes)
        {
            int i1 = 0;
            int i2 = 0;

            List<ReadRatio> list1 = ratios1.get(chromosome);
            List<ReadRatio> list2 = ratios2.get(chromosome);

            while (i1 < list1.size() && i2 < list2.size())
            {
                ReadRatio readRatio1 = list1.get(i1);
                ReadRatio readRatio2 = list2.get(i2);

                if(readRatio1.position() < readRatio2.position())
                {
                    mCombinedRatios.put(chromosome, readRatio1);
                    ++i1;
                }
                else if(readRatio1.position() > readRatio2.position())
                {
                    mCombinedRatios.put(chromosome, readRatio2);
                    ++i2;
                }
                else
                {
                    // shouldn't be here
                    CB_LOGGER.error("on target and off target ratio at same position: {}:{}",
                        readRatio1.chromosome(), readRatio1.position());

                    throw new RuntimeException("on target and off target ratio at same position");
                }
            }

            // put the rest in
            while(i1 < list1.size())
            {
                mCombinedRatios.put(chromosome, list1.get(i1));
                ++i1;
            }
            while(i2 < list2.size())
            {
                mCombinedRatios.put(chromosome, list2.get(i2));
                ++i2;
            }
        }
    }

    private void populateOnTargetRatios(
            final Map<GenomePosition, Double> targetRelativeEnrichment, final ListMultimap<Chromosome,ReadRatio> rawRatios)
    {
        // find all the ratios that are inside the target enriched regions
        // we filter out all the regions with 0 gc normalised ratios, as they do not actually
        // correctly reflect the amount of enrichment, and also very rare

        List<Double> targetRegionsGcRatios = Lists.newArrayList();

        for(Map.Entry<Chromosome,ReadRatio> entry : rawRatios.entries())
        {
            ReadRatio readRatio = entry.getValue();

            if(readRatio.ratio() <= 0)
                continue;

            if(!targetRelativeEnrichment.containsKey(readRatio))
                continue;

            targetRegionsGcRatios.add(readRatio.ratio());
        }

        double targetRegionGcRatioMedian = 1.0;

        if(!targetRegionsGcRatios.isEmpty())
        {
            targetRegionGcRatioMedian = Doubles.median(targetRegionsGcRatios);
        }
        else
        {
            CB_LOGGER.warn("target region gc ratios is empty");
        }

        CB_LOGGER.printf(Level.INFO, "targeted mode GC ratio median: %.3f", targetRegionGcRatioMedian);

        mOnTargetRatios.clear();

        for(Map.Entry<Chromosome,ReadRatio> entry : rawRatios.entries())
        {
            ReadRatio readRatio = entry.getValue();

            Double relativeEnrichment = targetRelativeEnrichment.get(readRatio);

            if(relativeEnrichment == null)
                continue;

            double enrichmentAdjRatio = readRatio.ratio() / targetRegionGcRatioMedian / relativeEnrichment;

            CB_LOGGER.debug(format("%s:%d relative enrichment: %.3f, on target ratio: %.3f",
                    readRatio.chromosome(), readRatio.position(), relativeEnrichment, enrichmentAdjRatio));

            mOnTargetRatios.put(entry.getKey(), ImmutableReadRatio.builder().from(readRatio).ratio(enrichmentAdjRatio).build());
        }
    }

    // we create a pan window ratio by taking the median count of super windows that combine multiple windows
    private void populateOffTargetRatios(
        final ListMultimap<Chromosome, ReadRatio> rawRatios, final int offTargetWindowSize, final List<GenomePosition> targetRegions)
    {
        ArrayListMultimap<Chromosome,ReadRatio> unnormalizedRatios = ArrayListMultimap.create();

        Window window = new Window(offTargetWindowSize);

        for(Chromosome chromosome : rawRatios.keySet())
        {
            int currentWindowStart = -1;

            List<Double> windowGcRatios = Lists.newArrayList();

            // we need this to make sure we get consistent chromosome name (1 vs chr1)
            String chromosomeStr = rawRatios.get(chromosome).get(0).chromosome();

            for(ReadRatio readRatio : rawRatios.get(chromosome))
            {
                // todo: make sure this is sorted

                int windowStart = window.start(readRatio.position());

                if(windowStart != currentWindowStart)
                {
                    if(currentWindowStart != -1)
                    {
                        ReadRatio unnormalizedRatio = unnormalizedOffTargetRatio(
                            offTargetWindowSize, chromosomeStr, currentWindowStart, windowGcRatios, targetRegions);

                        if(unnormalizedRatio != null)
                        {
                            unnormalizedRatios.put(chromosome, unnormalizedRatio);
                        }
                    }

                    currentWindowStart = windowStart;
                    windowGcRatios.clear();
                }

                if(readRatio.ratio() >= 0)
                    windowGcRatios.add(readRatio.ratio());
            }

            ReadRatio unnormalizedRatio = unnormalizedOffTargetRatio(
                offTargetWindowSize, chromosomeStr, currentWindowStart, windowGcRatios, targetRegions);

            if(unnormalizedRatio != null)
                unnormalizedRatios.put(chromosome, unnormalizedRatio);
        }

        // now we want to normalise all of those off target gc ratios by the median
        List<Double> values = Lists.newArrayList();
        unnormalizedRatios.values().forEach(x -> values.add(x.ratio()));
        double median = Doubles.median(values);

        CB_LOGGER.debug("normalizing {} off target windows ratio by median: {}", unnormalizedRatios.size(), median);

        mOffTargetRatios.clear();

        //for ((key, value) in unnormalizedRatios.entries())
        for(Map.Entry<Chromosome,ReadRatio> entry : unnormalizedRatios.entries())
        {
            ReadRatio readRatio = entry.getValue();
            double normalizedRatio = readRatio.ratio() / median;

            mOffTargetRatios.put(entry.getKey(), ImmutableReadRatio.builder().from(readRatio).ratio(normalizedRatio).build());
        }
    }

    @Nullable
    private static ReadRatio unnormalizedOffTargetRatio(
            int offTargetWindowSize, final String chromosome, int windowStart, final List<Double> windowGcRatios,
            final List<GenomePosition> targetRegions)
    {
        int minNumGcRatios = (int)round(offTargetWindowSize / CobaltConstants.WINDOW_SIZE * CobaltConstants.MIN_OFF_TARGET_WINDOW_RATIO);

        int windowEnd = windowStart + offTargetWindowSize - 1;

        // the window position is the middle
        int windowMid = windowStart + offTargetWindowSize / 2;

        if(windowGcRatios.size() < minNumGcRatios)
        {
            // if we don't have enough sub windows with valid values then we skip this
            CB_LOGGER.trace( "off target window: {}:{} ({} - {}), not enough sub window",
                chromosome, windowMid, windowStart, windowEnd);
            return null;
        }

        // check that this window does not contain any target regions
        for(GenomePosition targetRegion : targetRegions)
        {
            if(targetRegion.chromosome().equals(chromosome) && windowStart <= targetRegion.position() && windowEnd > targetRegion.position())
            {
                // this window contains a target region
                CB_LOGGER.trace("off target window: {}:{} ({} - {}), contains target region",
                    chromosome, windowMid, windowStart, windowEnd);
                return null;
            }
        }

        double median = Doubles.median(windowGcRatios);

        CB_LOGGER.debug("off target window: {}:{} ({} - {}), num sub windows: {}, median: {}",
                chromosome, windowMid, windowStart, windowEnd, windowGcRatios.size(), format("%.4f", median));

        return !Double.isNaN(median) ? ImmutableReadRatio.builder().chromosome(chromosome).position(windowMid).ratio(median).build() : null;
    }
}