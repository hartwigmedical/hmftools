package com.hartwig.hmftools.gripss.filters;

import static com.hartwig.hmftools.common.sv.LineElements.isMobileLineElement;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.fromRefAlt;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_BAQ;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_BQ;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.SV_FRAG_COUNT;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.QUAL;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsDouble;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.gripss.common.VcfUtils.sglFragmentCount;

import com.hartwig.hmftools.common.sv.VariantAltInsertCoords;
import com.hartwig.hmftools.common.variant.GenotypeIds;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class HardFilters
{
    private final FilterConstants mFilterConstants;
    private final boolean mGermlineMode;

    public HardFilters(final FilterConstants filterConstants, final boolean germlineMode)
    {
        mFilterConstants = filterConstants;
        mGermlineMode = germlineMode;
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
            final VariantAltInsertCoords altInsertCoords = fromRefAlt(variant.getAlleles().get(1).getDisplayString(), ref);

            boolean isLineInsertion = isMobileLineElement(altInsertCoords.Orient, altInsertCoords.InsertSequence);

            final String qualTag = isLineInsertion ? GRIDSS_BQ : GRIDSS_BAQ;
            qual = getGenotypeAttributeAsDouble(tumorGenotype, qualTag, 0);
        }
        else
        {
            qual = getGenotypeAttributeAsDouble(tumorGenotype, QUAL, 0);
        }

        return qual < mFilterConstants.MinTumorQual;
    }

    private boolean hasExcessiveReferenceSupport(final VariantContext variant, final GenotypeIds genotypeIds, boolean isSgl)
    {
        if(!genotypeIds.hasReference() || mGermlineMode)
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
            refFrags = getGenotypeAttributeAsInt(refGenotype, SV_FRAG_COUNT, 0);
            tumorFrags = getGenotypeAttributeAsInt(tumorGenotype, SV_FRAG_COUNT, 0);
        }

        if(refFrags > mFilterConstants.HardMaxNormalRelativeSupport * tumorFrags)
            return true;

        if(refFrags > mFilterConstants.HardMaxNormalAbsoluteSupport && refFrags > mFilterConstants.SoftMaxNormalRelativeSupport * tumorFrags)
            return true;

        return false;
    }


}
