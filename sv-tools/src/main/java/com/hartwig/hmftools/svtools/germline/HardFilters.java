package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.svtools.germline.CommonFilters.isPolyATSequence;
import static com.hartwig.hmftools.svtools.germline.FilterType.HARD_MIN_QUAL;
import static com.hartwig.hmftools.svtools.germline.FilterType.PASS;
import static com.hartwig.hmftools.svtools.germline.GermlineFilters.getDoubleValue;
import static com.hartwig.hmftools.svtools.germline.SvData.hasLength;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.QUAL;

import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;

import htsjdk.variant.variantcontext.VariantContext;

public class HardFilters
{
    private final FilterConstants mFilterConstants;
    private final HotspotCache mHotspotCache;

    public HardFilters(final FilterConstants filterConstants, final HotspotCache hotspotCache)
    {
        mFilterConstants = filterConstants;
        mHotspotCache = hotspotCache;
    }

    public boolean failsMinTumorQuality(final VariantContext variant, int tumorOrdinal)
    {
        double qual = getDoubleValue(variant.getGenotype(tumorOrdinal), QUAL);

        if(qual >= mFilterConstants.MinTumorQual)
            return false;

        if(StructuralVariantFactory.isSingleBreakend(variant))
            return true;

        // keep if matches a hotspot region (without knowing the other end yet)
        return mHotspotCache.matchesHotspotBreakend(variant.getContig(), variant.getStart());
    }

    public FilterType getFilterType(final StructuralVariant sv)
    {
        // the following hard filters are applied:
        // - no mate
        // - below min tumor qual
        // - excessive normal support

        FilterType filterType = PASS;

        if(belowMinQual(sv))
        {
            filterType = HARD_MIN_QUAL;
        }

        if(filterType != PASS)
        {
            // check hotspot rescue
            boolean canRescue = true;

            if(hasLength(sv.type()) && SvData.length(sv) < FilterConstants.SHORT_RESCUE_LENGTH)
                canRescue = false;
            else if(isPolyATSequence(sv.startContext()) || (sv.endContext() != null && isPolyATSequence(sv.endContext())))
                canRescue = false;

            if(canRescue && mHotspotCache.matchesHotspot(sv))
                filterType = PASS;
        }

        return filterType;
    }

    /*
        fun normalCoverageFilter(minNormalCoverage: Int): Boolean {
        if (normalGenotype == null) {
            return false
        }

        val supportingFragments = normalGenotype.fragmentSupport(isSingle)
        val ref = normalGenotype.refSupportRead()
        val refPair = normalGenotype.refSupportReadPair()

        return supportingFragments + ref + refPair < minNormalCoverage
    }

    fun normalSupportAbsoluteFilter(maxNormalAbsoluteSupport: Int): Boolean {
        if (normalGenotype == null) {
            return false
        }

        val normalSupport = normalGenotype.fragmentSupport(isSingle)
        return normalSupport > maxNormalAbsoluteSupport
    }

    fun normalSupportRelativeFilter(maxNormalRelativeSupport: Double): Boolean {
        if (normalGenotype == null ) {
            return false
        }

        val normalSupport = normalGenotype.fragmentSupport(isSingle)
        val tumorSupport = tumorGenotype.fragmentSupport(isSingle)

        return normalSupport > maxNormalRelativeSupport * tumorSupport
    }

     */

    private boolean belowMinQual(final StructuralVariant sv)
    {
        double qual = getDoubleValue(sv.startContext().getGenotype(0), QUAL);

        if(qual < mFilterConstants.MinTumorQual)
            return true;

        if(sv.endContext() != null)
        {
            qual = getDoubleValue(sv.endContext().getGenotype(0), QUAL);

            if(qual < mFilterConstants.MinTumorQual)
                return true;
        }

        return false;
    }

}
