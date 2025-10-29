package com.hartwig.hmftools.cobalt.targeted;

import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.cobalt.ratio.RatioMapper;

import org.apache.commons.lang3.Validate;

import tech.tablesaw.api.Table;

public class TargetedRegionsFilter implements RatioMapper
{
    private final Table mTargetRegionEnrichment;

    public TargetedRegionsFilter(final Table targetRegionEnrichment)
    {
        mTargetRegionEnrichment = targetRegionEnrichment;
    }

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
        return onTargetRatios;
    }


}