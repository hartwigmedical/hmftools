package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.svtools.germline.CommonFilters.isPolyATSequence;
import static com.hartwig.hmftools.svtools.germline.FilterConstants.MIN_QUAL_SUPPORT;
import static com.hartwig.hmftools.svtools.germline.FilterConstants.SHORT_RESCUE_LENGTH;
import static com.hartwig.hmftools.svtools.germline.FilterType.HARD_MIN_QUAL;
import static com.hartwig.hmftools.svtools.germline.FilterType.PASS;
import static com.hartwig.hmftools.svtools.germline.GermlineFilters.getDoubleValue;
import static com.hartwig.hmftools.svtools.germline.SvData.hasLength;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.QUAL;

import com.hartwig.hmftools.common.sv.StructuralVariant;

public class HardFilters
{
    private final HotspotCache mHotspotCache;

    public HardFilters(final HotspotCache hotspotCache)
    {
        mHotspotCache = hotspotCache;
    }

    public FilterType getFilterType(final StructuralVariant sv)
    {
        FilterType filterType = PASS;

        if(belowMinQual(sv))
        {
            filterType = HARD_MIN_QUAL;
        }

        if(filterType != PASS)
        {
            // check hotspot rescue
            boolean canRescue = true;

            if(hasLength(sv.type()) && SvData.length(sv) < SHORT_RESCUE_LENGTH)
                canRescue = false;
            else if(isPolyATSequence(sv.startContext()) || (sv.endContext() != null && isPolyATSequence(sv.endContext())))
                canRescue = false;

            if(canRescue && mHotspotCache.matchesHotspot(sv))
                filterType = PASS;
        }

        return filterType;
    }

    private static boolean belowMinQual(final StructuralVariant sv)
    {
        double qual = getDoubleValue(sv.startContext().getGenotype(0), QUAL);

        if(qual < MIN_QUAL_SUPPORT)
            return true;

        if(sv.endContext() != null)
        {
            qual = getDoubleValue(sv.endContext().getGenotype(0), QUAL);

            if(qual < MIN_QUAL_SUPPORT)
                return true;
        }

        return false;
    }

}
