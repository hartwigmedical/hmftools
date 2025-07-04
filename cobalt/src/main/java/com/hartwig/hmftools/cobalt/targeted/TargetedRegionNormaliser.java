package com.hartwig.hmftools.cobalt.targeted;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;

import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.cobalt.ratio.RatioMapper;

import org.apache.logging.log4j.Level;

import tech.tablesaw.api.DoubleColumn;
import tech.tablesaw.api.Table;
import tech.tablesaw.columns.numbers.NumberPredicates;

public class TargetedRegionNormaliser implements RatioMapper
{
    @Override
    public Table mapRatios(final Table onTargetRatios)
    {
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
}