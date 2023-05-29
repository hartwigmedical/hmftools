package com.hartwig.hmftools.cobalt.targeted;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltUtils.replaceColumn;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.cobalt.ChromosomePositionCodec;
import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.cobalt.CobaltConstants;
import com.hartwig.hmftools.cobalt.lowcov.LowCoverageRatioBuilder;
import com.hartwig.hmftools.cobalt.ratio.GcNormalizedRatioBuilder;
import com.hartwig.hmftools.cobalt.ratio.RatioBuilder;
import com.hartwig.hmftools.common.cobalt.ImmutableReadRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.lang3.Validate;
import org.jetbrains.annotations.Nullable;

import tech.tablesaw.api.Table;
import tech.tablesaw.columns.numbers.NumberPredicates;

public class TargetedRatioBuilder implements RatioBuilder
{
    private Table mOnTargetRatios;
    private Table mOffTargetRatios;
    private Table mCombinedRatios;

    private final ChromosomePositionCodec mChromosomePosCodec;

    public TargetedRatioBuilder(final Table rawRatios, final Table targetRegionEnrichment,
            ChromosomePositionCodec chromosomePosCodec)
    {
        mOnTargetRatios = Table.create();
        mOffTargetRatios = Table.create();
        mCombinedRatios = Table.create();
        mChromosomePosCodec = chromosomePosCodec;

        Validate.isTrue(rawRatios.longColumn(CobaltColumns.ENCODED_CHROMOSOME_POS).isMissing().isEmpty());
        Validate.isTrue(targetRegionEnrichment.longColumn(CobaltColumns.ENCODED_CHROMOSOME_POS).isMissing().isEmpty());

        // merge in the targeted region columns
        Table ratiosWithTargetRegion = rawRatios.joinOn(CobaltColumns.ENCODED_CHROMOSOME_POS).inner(targetRegionEnrichment);

        // resort it, the join messes up with the ordering
        ratiosWithTargetRegion = ratiosWithTargetRegion.sortAscendingOn(CobaltColumns.ENCODED_CHROMOSOME_POS);

        CB_LOGGER.info("ratios with targetd region: {}", ratiosWithTargetRegion);

        populateOnTargetRatios(ratiosWithTargetRegion, targetRegionEnrichment);
        populateOffTargetRatios(ratiosWithTargetRegion, CobaltConstants.OFF_TARGET_WINDOW_SIZE, targetRegionEnrichment);
        CB_LOGGER.info("{} on target GC ratios, {} off target GC ratios", mOnTargetRatios.rowCount(), mOffTargetRatios.rowCount());

        // populateCombinedRatios(mOnTargetRatios, mOffTargetRatios);
    }

    // we use on target ratios only for now
    @Override
    public Table ratios() { return mOffTargetRatios; }

    public Table onTargetRatios() { return mOnTargetRatios; }
    public Table offTargetRatios() { return mOffTargetRatios; }

    private void populateCombinedRatios(final Table ratios1, final Table ratios2)
    {
        mCombinedRatios = ratios1.append(ratios2);

        /*Set<Chromosome> chromosomes = Sets.newIdentityHashSet();

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
         */
    }

    private void populateOnTargetRatios(final Table rawRatios, final Table targetRelativeEnrichment)
    {
        // find all the ratios that are inside the target enriched regions
        // we filter out all the regions with 0 gc normalised ratios, as they do not actually
        // correctly reflect the amount of enrichment, and also very rare

        /*
        List<Double> targetRegionsGcRatios = new ArrayList<>();

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
         */

        Table onTargetRatios = rawRatios.copy();
        double median = onTargetRatios.doubleColumn("ratio").filter(NumberPredicates.isNonNegative).median();

        // normalise the ratio by relative enrichment
        replaceColumn(onTargetRatios, "ratio",
                onTargetRatios.doubleColumn("ratio")
                        .divide(onTargetRatios.doubleColumn("relativeEnrichment")));

        // apply gc normalisation
        GcNormalizedRatioBuilder gcNormalizedRatioBuilder = new GcNormalizedRatioBuilder(onTargetRatios, false);
        mOnTargetRatios = gcNormalizedRatioBuilder.ratios();

        /*
        mOnTargetRatios = rawRatios.copy();

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
         */
    }

    // we create a pan window ratio by taking the median count of super windows that combine multiple windows
    private void populateOffTargetRatios(final Table ratiosWithTargetRegion, final int offTargetWindowSize,
            final Table targetRegions)
    {
        Table offTargetRatios = ratiosWithTargetRegion.where(
                ratiosWithTargetRegion.booleanColumn("offTarget").asSelection()
                        .and(ratiosWithTargetRegion.doubleColumn("ratio").isNonNegative())
                        .and(ratiosWithTargetRegion.doubleColumn("relativeEnrichment").isNotMissing()));

        // double median = offTargetRatios.doubleColumn("ratio").median();

        // normalise the ratio by relative enrichment
        replaceColumn(offTargetRatios, "ratio",
                offTargetRatios.doubleColumn("ratio")
                        .divide(offTargetRatios.doubleColumn("relativeEnrichment")));

        CB_LOGGER.info("off target after enrichment normalisation: \n{}", offTargetRatios);

        // apply gc normalisation
        GcNormalizedRatioBuilder gcNormalizedRatioBuilder = new GcNormalizedRatioBuilder(offTargetRatios, true);
        offTargetRatios = gcNormalizedRatioBuilder.ratios();

        CB_LOGGER.info("off target gc normalisation: \n{}", gcNormalizedRatioBuilder.gcMedianReadCountTable());
        CB_LOGGER.info("off target after gc normalisation: \n{}", offTargetRatios);

        // next we do low coverage
        offTargetRatios = new LowCoverageRatioBuilder(offTargetRatios, 500, mChromosomePosCodec).ratios();

        CB_LOGGER.info("off target after consolidation: \n{}", offTargetRatios);

        mOffTargetRatios = offTargetRatios;

        // remove any with invalid ratios
        // mOffTargetRatios = offTargetRatios.where(offTargetRatios.doubleColumn(CobaltColumns.RATIO).)

        /*
        Window window = new Window(offTargetWindowSize);

        for(String chromosome : rawRatios.stringColumn("chromosome").unique())
        {
            int currentWindowStart = 1;

            List<ReadRatio> windowGcRatios = new ArrayList<>();

            // we need this to make sure we get consistent chromosome name (1 vs chr1)


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
                    windowGcRatios.add(readRatio);
            }

            ReadRatio unnormalizedRatio = unnormalizedOffTargetRatio(
                offTargetWindowSize, chromosomeStr, currentWindowStart, windowGcRatios, targetRegions);

            if(unnormalizedRatio != null)
                unnormalizedRatios.put(chromosome, unnormalizedRatio);
        }

        // now we want to normalise all of those off target gc ratios by the median
        List<Double> values = new ArrayList<>();
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
         */
    }

    @Nullable
    private static ReadRatio unnormalizedOffTargetRatio(
            int offTargetWindowSize, final String chromosome, int windowStart, final List<ReadRatio> windowGcRatios,
            final List<GenomePosition> targetRegions)
    {
        int minNumGcRatios = (int)round(offTargetWindowSize / CobaltConstants.WINDOW_SIZE * CobaltConstants.MIN_OFF_TARGET_WINDOW_RATIO);

        int windowEnd = windowStart + offTargetWindowSize - 1;

        // the window position is the middle
        int windowMid = windowStart + offTargetWindowSize / 2;

        // check for targeted regions, we want to remove them
        for(GenomePosition targetRegion : targetRegions)
        {
            if(targetRegion.chromosome().equals(chromosome) && windowStart <= targetRegion.position() && windowEnd > targetRegion.position())
            {
                // this window contains a target region
                int removeStart = targetRegion.position() - 2 * CobaltConstants.WINDOW_SIZE;
                int removeEnd = targetRegion.position() + 2 * CobaltConstants.WINDOW_SIZE;
                CB_LOGGER.trace("off target window: {}:{} ({} - {}), contains target region",
                    chromosome, windowMid, windowStart, windowEnd);

                windowGcRatios.removeIf(o -> o.position() >= removeStart && o.position() <= removeEnd);
            }
        }

        if(windowGcRatios.size() < minNumGcRatios)
        {
            // if we don't have enough sub windows with valid values then we skip this
            CB_LOGGER.trace( "off target window: {}:{} ({} - {}), not enough sub window",
                    chromosome, windowMid, windowStart, windowEnd);
            return null;
        }

        double median = Doubles.median(windowGcRatios.stream().map(ReadRatio::ratio).collect(Collectors.toList()));

        CB_LOGGER.debug("off target window: {}:{} ({} - {}), num sub windows: {}, median: {}",
                chromosome, windowMid, windowStart, windowEnd, windowGcRatios.size(), format("%.4f", median));

        return !Double.isNaN(median) ? ImmutableReadRatio.builder().chromosome(chromosome).position(windowMid).ratio(median).build() : null;
    }
}