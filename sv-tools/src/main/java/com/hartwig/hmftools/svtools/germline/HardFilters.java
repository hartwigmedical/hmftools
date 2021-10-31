package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.svtools.germline.CommonFilters.isPolyATSequence;
import static com.hartwig.hmftools.svtools.germline.GermlineFilters.getDoubleValue;
import static com.hartwig.hmftools.svtools.germline.SvData.hasLength;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.BVF;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.QUAL;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.VF;

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
        double qual = getDoubleValue(variant.getGenotype(genotypeIds.TumorOrdinal), QUAL);
        return qual < mFilterConstants.MinTumorQual;
    }

    private boolean hasExcessiveReferenceSupport(final VariantContext variant, final GenotypeIds genotypeIds)
    {
        if(!genotypeIds.hasReference())
            return false;

        final String supportTag = StructuralVariantFactory.isSingleBreakend(variant) ? BVF : VF;
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
