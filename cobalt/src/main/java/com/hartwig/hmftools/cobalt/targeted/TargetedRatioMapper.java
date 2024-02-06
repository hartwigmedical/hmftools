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
import com.hartwig.hmftools.cobalt.lowcov.LowCoverageRatioMapper;
import com.hartwig.hmftools.cobalt.ratio.GcNormalizedRatioMapper;
import com.hartwig.hmftools.cobalt.ratio.RatioMapper;
import com.hartwig.hmftools.common.cobalt.ImmutableReadRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.lang3.Validate;
import org.apache.logging.log4j.Level;
import org.jetbrains.annotations.Nullable;

import tech.tablesaw.api.DoubleColumn;
import tech.tablesaw.api.Table;
import tech.tablesaw.columns.numbers.NumberPredicates;

public class TargetedRatioMapper implements RatioMapper
{
    private final Table mTargetRegionEnrichment;
    private final ChromosomePositionCodec mChromosomePosCodec;

    public TargetedRatioMapper(final Table targetRegionEnrichment,
            ChromosomePositionCodec chromosomePosCodec)
    {
        mTargetRegionEnrichment = targetRegionEnrichment;
        mChromosomePosCodec = chromosomePosCodec;
    }

    // we use on target ratios only for now
    @Override
    public Table mapRatios(final Table inputRatios)
    {
        return onTargetRatios(inputRatios);
    }

    private void populateCombinedRatios(final Table ratios1, final Table ratios2)
    {
        // mCombinedRatios = ratios1.append(ratios2);
    }

    Table onTargetRatios(final Table inputRatios)
    {
        // find all the ratios that are inside the target enriched regions
        // we filter out all the regions with 0 gc normalised ratios, as they do not actually
        // correctly reflect the amount of enrichment, and also very rare

        Validate.isTrue(inputRatios.longColumn(CobaltColumns.ENCODED_CHROMOSOME_POS).isMissing().isEmpty());
        Validate.isTrue(mTargetRegionEnrichment.longColumn(CobaltColumns.ENCODED_CHROMOSOME_POS).isMissing().isEmpty());

        // merge in the targeted region columns
        Table onTargetRatios = inputRatios.joinOn(CobaltColumns.ENCODED_CHROMOSOME_POS).inner(mTargetRegionEnrichment);

        // resort it, the join messes up with the ordering
        onTargetRatios = onTargetRatios.sortAscendingOn(CobaltColumns.ENCODED_CHROMOSOME_POS);

        double targetRegionGcRatioMedian = onTargetRatios.doubleColumn(CobaltColumns.RATIO).filter(NumberPredicates.isNonNegative).median();

        CB_LOGGER.printf(Level.INFO, "targeted mode GC ratio median: %.3f", targetRegionGcRatioMedian);

        // normalise the ratio by relative enrichment and targeted region median
        DoubleColumn onTargetRatioColumn = onTargetRatios.doubleColumn(CobaltColumns.RATIO)
                .divide(onTargetRatios.doubleColumn("relativeEnrichment"))
                .divide(targetRegionGcRatioMedian)
                .map(d -> Double.isFinite(d) ? d : Double.NaN) // protect against division by 0
                .setName(CobaltColumns.RATIO);

        onTargetRatios.replaceColumn(onTargetRatioColumn);

        return onTargetRatios;
    }

    // we create a pan window ratio by taking the median count of super windows that combine multiple windows
    Table offTargetRatios(final Table inputRatios)
    {
        // merge in the targeted region columns
        Table offTargetRatios = inputRatios.joinOn(CobaltColumns.ENCODED_CHROMOSOME_POS).inner(mTargetRegionEnrichment);

        // resort it, the join messes up with the ordering
        offTargetRatios = offTargetRatios.sortAscendingOn(CobaltColumns.ENCODED_CHROMOSOME_POS);

        offTargetRatios = offTargetRatios.where(
                        offTargetRatios.booleanColumn("offTarget").asSelection()
                        .and(offTargetRatios.doubleColumn("ratio").isNonNegative())
                        .and(offTargetRatios.doubleColumn("relativeEnrichment").isNotMissing()));

        // double median = offTargetRatios.doubleColumn("ratio").median();

        // normalise the ratio by relative enrichment
        replaceColumn(offTargetRatios, "ratio",
                offTargetRatios.doubleColumn("ratio")
                        .divide(offTargetRatios.doubleColumn("relativeEnrichment")));

        CB_LOGGER.info("off target after enrichment normalisation: \n{}", offTargetRatios);

        // next we do low coverage
        offTargetRatios = new LowCoverageRatioMapper(1000, mChromosomePosCodec).mapRatios(offTargetRatios);

        CB_LOGGER.info("off target after consolidation: \n{}", offTargetRatios);

        // apply gc normalisation
        GcNormalizedRatioMapper gcNormalizedRatioMapper = new GcNormalizedRatioMapper();
        offTargetRatios = gcNormalizedRatioMapper.mapRatios(offTargetRatios);

        CB_LOGGER.info("off target gc normalisation: \n{}", gcNormalizedRatioMapper.gcMedianReadDepthTable());
        CB_LOGGER.info("off target after gc normalisation: \n{}", offTargetRatios);

        return offTargetRatios;

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
        int minNumGcRatios = (int)round((double)offTargetWindowSize / CobaltConstants.WINDOW_SIZE * CobaltConstants.MIN_OFF_TARGET_WINDOW_RATIO);

        int windowEnd = windowStart + offTargetWindowSize - 1;

        // the window position is the middle
        int windowMid = windowStart + offTargetWindowSize / 2;

        // check for targeted regions, we want to remove them
        for (GenomePosition targetRegion : targetRegions)
        {
            if (targetRegion.chromosome().equals(chromosome) && windowStart <= targetRegion.position() && windowEnd > targetRegion.position())
            {
                // this window contains a target region
                int removeStart = targetRegion.position() - 2 * CobaltConstants.WINDOW_SIZE;
                int removeEnd = targetRegion.position() + 2 * CobaltConstants.WINDOW_SIZE;
                CB_LOGGER.trace("off target window: {}:{} ({} - {}), contains target region",
                    chromosome, windowMid, windowStart, windowEnd);

                windowGcRatios.removeIf(o -> o.position() >= removeStart && o.position() <= removeEnd);
            }
        }

        if (windowGcRatios.size() < minNumGcRatios)
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