package com.hartwig.hmftools.gripss.filters;

import static com.hartwig.hmftools.common.sv.LineElements.isMobileLineElement;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsDouble;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.gripss.common.VariantAltInsertCoords.parseRefAlt;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_BAQ;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_BQ;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_QUAL;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_VF;
import static com.hartwig.hmftools.gripss.common.VcfUtils.sglFragmentCount;

import com.hartwig.hmftools.gripss.common.VariantAltInsertCoords;
import com.hartwig.hmftools.gripss.common.GenotypeIds;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class HardFilters
{
    private final FilterConstants mFilterConstants;

    public HardFilters(final FilterConstants filterConstants)
    {
        mFilterConstants = filterConstants;
    }

    public boolean isFiltered(final VariantContext variant, final GenotypeIds genotypeIds, boolean isSgl)
    {
        if(isSgl && mFilterConstants.FilterSGLs)
            return true;

        // the following hard filters are applied:
        // - below min tumor qual
        // - excessive normal support
        if(belowMinQual(variant, genotypeIds, isSgl))
            return true;

        if(hasExcessiveReferenceSupport(variant, genotypeIds, isSgl))
            return true;

        return false;
    }

    public boolean belowMinQual(final VariantContext variant, final GenotypeIds genotypeIds, boolean isSgl)
    {
        Genotype tumorGenotype = variant.getGenotype(genotypeIds.TumorOrdinal);
        double qual;

        if(isSgl)
        {
            String ref = variant.getAlleles().get(0).getDisplayString();
            final VariantAltInsertCoords altInsertCoords = parseRefAlt(variant.getAlleles().get(1).getDisplayString(), ref);

            boolean isLineInsertion = isMobileLineElement(altInsertCoords.Orientation, altInsertCoords.InsertSequence);

            final String qualTag = isLineInsertion ? VT_BQ : VT_BAQ;
            qual = getGenotypeAttributeAsDouble(tumorGenotype, qualTag, 0);
        }
        else
        {
            qual = getGenotypeAttributeAsDouble(tumorGenotype, VT_QUAL, 0);
        }

        return qual < mFilterConstants.MinTumorQual;
    }

    private boolean hasExcessiveReferenceSupport(final VariantContext variant, final GenotypeIds genotypeIds, boolean isSgl)
    {
        if(!genotypeIds.hasReference() || genotypeIds.GermlineMode)
            return false;

        Genotype refGenotype = variant.getGenotype(genotypeIds.ReferenceOrdinal);
        Genotype tumorGenotype = variant.getGenotype(genotypeIds.TumorOrdinal);

        int refFrags;
        int tumorFrags;

        if(isSgl)
        {
            refFrags = sglFragmentCount(refGenotype);
            tumorFrags = sglFragmentCount(tumorGenotype);
        }
        else
        {
            refFrags = getGenotypeAttributeAsInt(refGenotype, VT_VF, 0);
            tumorFrags = getGenotypeAttributeAsInt(tumorGenotype, VT_VF, 0);
        }

        if(refFrags > mFilterConstants.HardMaxNormalRelativeSupport * tumorFrags)
            return true;

        if(refFrags > mFilterConstants.HardMaxNormalAbsoluteSupport && refFrags > mFilterConstants.SoftMaxNormalRelativeSupport * tumorFrags)
            return true;

        return false;
    }


}
