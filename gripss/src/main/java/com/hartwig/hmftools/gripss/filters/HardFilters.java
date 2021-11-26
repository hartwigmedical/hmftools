package com.hartwig.hmftools.gripss.filters;

import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_BVF;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_QUAL;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_VF;

import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.gripss.common.VcfUtils;
import com.hartwig.hmftools.gripss.common.GenotypeIds;

import htsjdk.variant.variantcontext.VariantContext;

public class HardFilters
{
    private final FilterConstants mFilterConstants;

    public HardFilters(final FilterConstants filterConstants)
    {
        mFilterConstants = filterConstants;
    }

    public boolean isFiltered(final VariantContext variant, final GenotypeIds genotypeIds)
    {
        // the following hard filters are applied:
        // - below min tumor qual
        // - excessive normal support

        if(belowMinQual(variant, genotypeIds))
            return true;

        if(hasExcessiveReferenceSupport(variant, genotypeIds))
            return true;

        return false;
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


}
