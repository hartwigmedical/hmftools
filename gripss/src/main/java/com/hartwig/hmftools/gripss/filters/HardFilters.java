package com.hartwig.hmftools.gripss.filters;

import static com.hartwig.hmftools.gripss.filters.CommonFilters.isPolyATSequence;
import static com.hartwig.hmftools.gripss.SvData.hasLength;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_BVF;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_QUAL;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_VF;

import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.gripss.GenotypeIds;
import com.hartwig.hmftools.gripss.HotspotCache;
import com.hartwig.hmftools.gripss.SvData;
import com.hartwig.hmftools.gripss.VcfUtils;

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

    public boolean isFiltered(final VariantContext variant, final GenotypeIds genotypeIds)
    {
        // the following hard filters are applied:
        // - below min tumor qual
        // - excessive normal support

        boolean failsFilters = false;

        if(belowMinQual(variant, genotypeIds))
        {
            failsFilters = true;
        }
        else if(hasExcessiveReferenceSupport(variant, genotypeIds))
        {
            failsFilters = true;
        }

        if(!failsFilters)
            return false;

        if(StructuralVariantFactory.isSingleBreakend(variant))
            return true;

        // keep if matches a hotspot region (without knowing the other end yet)
        return !mHotspotCache.matchesHotspotBreakend(variant.getContig(), variant.getStart());
    }

    private boolean belowMinQual(final VariantContext variant, final GenotypeIds genotypeIds)
    {
        double qual = VcfUtils.getGenotypeAttributeAsDouble(variant.getGenotype(genotypeIds.TumorOrdinal), VT_QUAL, 0);
        return qual < mFilterConstants.MinTumorQual;
    }

    private boolean hasExcessiveReferenceSupport(final VariantContext variant, final GenotypeIds genotypeIds)
    {
        if(!genotypeIds.hasReference())
            return false;

        final String supportTag = StructuralVariantFactory.isSingleBreakend(variant) ? VT_BVF : VT_VF;
        int refFrags = VcfUtils.getGenotypeAttributeAsInt(variant.getGenotype(genotypeIds.ReferenceOrdinal), supportTag, 0);
        int tumorFrags = VcfUtils.getGenotypeAttributeAsInt(variant.getGenotype(genotypeIds.TumorOrdinal), supportTag, 0);

        if(refFrags > mFilterConstants.HardMaxNormalRelativeSupport * tumorFrags)
            return true;

        if(refFrags > mFilterConstants.HardMaxNormalAbsoluteSupport && refFrags > mFilterConstants.SoftMaxNormalRelativeSupport * tumorFrags)
            return true;

        return false;
    }

    public boolean keepHotspotVariant(final StructuralVariant sv)
    {
        // check hotspot rescue
        if(hasLength(sv.type()) && SvData.length(sv) < FilterConstants.SHORT_RESCUE_LENGTH)
            return false;
        else if(isPolyATSequence(sv.startContext()) || (sv.endContext() != null && isPolyATSequence(sv.endContext())))
            return false;

        return mHotspotCache.matchesHotspot(sv);
    }
}
