package com.hartwig.hmftools.cobalt.targeted;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;

import com.hartwig.hmftools.cobalt.ChromosomePositionCodec;
import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.cobalt.ratio.RatioMapper;

import org.apache.commons.lang3.Validate;
import org.apache.logging.log4j.Level;

import tech.tablesaw.api.DoubleColumn;
import tech.tablesaw.api.Table;
import tech.tablesaw.columns.numbers.NumberPredicates;

public class TargetedRatioMapper implements RatioMapper
{
    private final Table mTargetRegionEnrichment;

    public TargetedRatioMapper(final Table targetRegionEnrichment)
    {
        mTargetRegionEnrichment = targetRegionEnrichment;
    }

    // we use on target ratios only for now
    @Override
    public Table mapRatios(final Table inputRatios)
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
}