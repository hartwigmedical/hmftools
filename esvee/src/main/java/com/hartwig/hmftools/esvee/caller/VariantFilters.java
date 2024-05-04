package com.hartwig.hmftools.esvee.caller;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sv.LineElements.POLY_A_HOMOLOGY;
import static com.hartwig.hmftools.common.sv.LineElements.POLY_T_HOMOLOGY;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOMSEQ;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.QUAL;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SPLIT_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.STRAND_BIAS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsDouble;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.HOM_INV_LENGTH;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.SGL_INS_SEQ_MIN_LENGTH;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.SGL_MAX_STRAND_BIAS;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.SGL_MIN_STRAND_BIAS;
import static com.hartwig.hmftools.esvee.common.FilterType.MAX_HOM_LENGTH_SHORT_INV;
import static com.hartwig.hmftools.esvee.common.FilterType.MAX_NORMAL_RELATIVE_SUPPORT;
import static com.hartwig.hmftools.esvee.common.FilterType.MAX_POLY_A_HOM_LENGTH;
import static com.hartwig.hmftools.esvee.common.FilterType.MIN_LENGTH;
import static com.hartwig.hmftools.esvee.common.FilterType.MIN_NORMAL_COVERAGE;
import static com.hartwig.hmftools.esvee.common.FilterType.MIN_QUALITY;
import static com.hartwig.hmftools.esvee.common.FilterType.MIN_TUMOR_AF;
import static com.hartwig.hmftools.esvee.common.FilterType.MODIFIED_AF;
import static com.hartwig.hmftools.esvee.common.FilterType.SGL_INSERT_SEQ_MIN_LENGTH;
import static com.hartwig.hmftools.esvee.common.FilterType.SGL_STRAND_BIAS;
import static com.hartwig.hmftools.esvee.common.FilterType.SHORT_SR_NORMAL;
import static com.hartwig.hmftools.esvee.common.FilterType.SHORT_SR_SUPPORT;
import static com.hartwig.hmftools.esvee.common.FilterType.SHORT_STRAND_BIAS;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.esvee.common.FilterType;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class VariantFilters
{
    private final FilterConstants mFilterConstants;
    private final boolean mGermlineMode;

    public VariantFilters(final FilterConstants filterConstants, final boolean germlineMode)
    {
        mFilterConstants = filterConstants;
        mGermlineMode = germlineMode;
    }

    public void applyFilters(final SvData sv, final FilterCache filterCache)
    {
        if(filterCache.isHotspot(sv))
        {
            // only limited filters are checked for hotspots
            if(modifiedAF(sv, true))
            {
                filterCache.addBreakendFilters(sv.breakendStart(), Lists.newArrayList(MODIFIED_AF));

                if(!sv.isSgl())
                    filterCache.addBreakendFilters(sv.breakendEnd(), Lists.newArrayList(MODIFIED_AF));
            }

            return;
        }

        List<FilterType> beStartFilters = null;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(sv.isSgl() && se == SE_END)
                continue;

            List<FilterType> filters = Lists.newArrayList();

            if(se == SE_START)
                beStartFilters = filters;

            Breakend breakend = sv.breakends()[se];

            if(normalCoverage(breakend))
                filters.add(MIN_NORMAL_COVERAGE);

            if(normalRelativeSupport(breakend))
                filters.add(MAX_NORMAL_RELATIVE_SUPPORT);

            if(allelicFrequency(sv, breakend))
                filters.add(MIN_TUMOR_AF);

            if(minQuality(sv, breakend))
                filters.add(MIN_QUALITY);

            if(shortSplitReadTumor(sv, breakend))
                filters.add(SHORT_SR_SUPPORT);

            if(shortSplitReadNormal(sv, breakend))
                filters.add(SHORT_SR_NORMAL);

            if(singleStrandBias(breakend))
                filters.add(SGL_STRAND_BIAS);

            if(singleInsertSequenceMinLength(breakend))
                filters.add(SGL_INSERT_SEQ_MIN_LENGTH);

            if(strandBias(sv, breakend))
                filters.add(SHORT_STRAND_BIAS);

            if((se == SE_END && beStartFilters.contains(MAX_POLY_A_HOM_LENGTH)) || polyATHomology(sv))
                filters.add(MAX_POLY_A_HOM_LENGTH);

            if((se == SE_END && beStartFilters.contains(MAX_HOM_LENGTH_SHORT_INV)) || homologyLengthFilterShortInversion(sv))
                filters.add(MAX_HOM_LENGTH_SHORT_INV);

            if((se == SE_END && beStartFilters.contains(MIN_LENGTH)) || minLength(sv))
                filters.add(MIN_LENGTH);

            if((se == SE_END && beStartFilters.contains(MODIFIED_AF)) || modifiedAF(sv, false))
                filters.add(MODIFIED_AF);

            if(!filters.isEmpty())
                filterCache.addBreakendFilters(breakend, filters);
        }
    }

    private boolean normalCoverage(final Breakend breakend)
    {
        if(breakend.RefGenotype == null || mGermlineMode)
            return false;

        int refSupportReads = getGenotypeAttributeAsInt(breakend.RefGenotype, REF_DEPTH, 0);
        int refSupportReadPairs = getGenotypeAttributeAsInt(breakend.RefGenotype, REF_DEPTH_PAIR, 0);

        return breakend.ReferenceFragments + refSupportReads + refSupportReadPairs < mFilterConstants.MinNormalCoverage;
    }

    private boolean normalRelativeSupport(final Breakend breakend)
    {
        if(breakend.RefGenotype == null || mGermlineMode)
            return false;

        return breakend.ReferenceFragments > mFilterConstants.SoftMaxNormalRelativeSupport * breakend.TumorFragments;
    }

    private boolean allelicFrequency(final SvData sv, final Breakend breakend)
    {
        double afThreshold = sv.isSgl() ? mFilterConstants.MinTumorAfBreakend : mFilterConstants.MinTumorAfBreakpoint;
        return breakend.allelicFrequency() < afThreshold;
    }

    private boolean modifiedAF(final SvData sv, boolean isHotspot)
    {
        double afLimit = isHotspot ? mFilterConstants.ModifiedAfHotspot : mFilterConstants.ModifiedAf;

        if(afLimit == 0)
            return false;

        double modifiedAf = calcModifiedAf(sv.breakendStart());

        if(!sv.isSgl())
        {
            modifiedAf = min(modifiedAf, calcModifiedAf(sv.breakendEnd()));
        }

        return modifiedAf < afLimit;
    }

    private static double calcModifiedAf(final Breakend breakend)
    {
        int readPairSupport = (breakend.isSgl() || !breakend.sv().isShortLocal()) ?
                getGenotypeAttributeAsInt(breakend.TumorGenotype, REF_DEPTH_PAIR, 0) : 0;

        int refSupport = getGenotypeAttributeAsInt(breakend.TumorGenotype, REF_DEPTH, 0);
        double totalSupport = breakend.TumorFragments + refSupport + readPairSupport;

        return totalSupport > 0 ? breakend.TumorFragments / totalSupport : 0;
    }

    private boolean minQuality(final SvData sv, final Breakend breakend)
    {
        double qualThreshold = sv.isSgl() ? mFilterConstants.MinQualBreakend : mFilterConstants.MinQualBreakpoint;

        if(mFilterConstants.LowQualRegion.containsPosition(breakend.Chromosome, breakend.Position))
            qualThreshold *= 0.5;

        return breakend.Qual < qualThreshold;
    }

    private boolean polyATHomology(final SvData sv)
    {
        String homologySequence = sv.contextStart().getAttributeAsString(HOMSEQ, "");
        return homologySequence.contains(POLY_A_HOMOLOGY) || homologySequence.contains(POLY_T_HOMOLOGY);
    }

    private static double calcStrandBias(final VariantContext variantContext)
    {
        double strandBias = variantContext.getAttributeAsDouble(STRAND_BIAS, 0.5);
        return max(strandBias, 1 - strandBias);
    }

    private boolean singleStrandBias(final Breakend breakend)
    {
        if(!breakend.isSgl() || breakend.IsLineInsertion)
            return false;

        if(mFilterConstants.LowQualRegion.containsPosition(breakend.Chromosome, breakend.Position))
            return false;

        double strandBias = calcStrandBias(breakend.Context);
        return strandBias < SGL_MIN_STRAND_BIAS || strandBias > SGL_MAX_STRAND_BIAS;
    }

    private boolean singleInsertSequenceMinLength(final Breakend breakend)
    {
        if(!breakend.isSgl() || breakend.IsLineInsertion)
            return false;

        return breakend.InsertSequence.length() < SGL_INS_SEQ_MIN_LENGTH;
    }

    private boolean homologyLengthFilterShortInversion(final SvData sv)
    {
        if(sv.type() != INV)
            return false;

        return sv.length() <= HOM_INV_LENGTH && sv.startHomology().length() > mFilterConstants.MaxHomLengthShortInv;
    }

    private boolean shortSplitReadTumor(final SvData sv, final Breakend breakend)
    {
        return sv.isShortLocal() && getSplitReadCount(breakend.TumorGenotype) == 0;
    }

    private boolean shortSplitReadNormal(final SvData sv, final Breakend breakend)
    {
        if(breakend.RefGenotype == null || mGermlineMode)
            return false;

        return sv.isShortLocal() && getSplitReadCount(breakend.RefGenotype) > 0;
    }

    private static int getSplitReadCount(final Genotype genotype)
    {
        return getGenotypeAttributeAsInt(genotype, SPLIT_FRAGS, 0);
    }

    private boolean strandBias(final SvData sv, final Breakend breakend)
    {
        if(sv.isShortLocal())
        {
            double strandBias = breakend.Context.getAttributeAsDouble(STRAND_BIAS, 0.5);
            return max(strandBias, 1 - strandBias) > mFilterConstants.MaxShortStrandBias;
        }

        return false;
    }

    private boolean minLength(final SvData sv)
    {
        if(sv.type() == DEL)
            return sv.length() + sv.insertSequence().length() - 1 < mFilterConstants.MinLength;
        else if(sv.type() == DUP)
            return sv.length() + sv.insertSequence().length() < mFilterConstants.MinLength;
        else if(sv.type() == INS)
            return sv.length() + sv.insertSequence().length() + 1 < mFilterConstants.MinLength;
        else
            return false;
    }


    // OLD HARD-FILTERS
    public boolean isHardFiltered(final VariantContext variant, final GenotypeIds genotypeIds, boolean isSgl)
    {
        if(isSgl && mFilterConstants.FilterSGLs)
            return true;

        // the following hard filters are applied:
        // - below min tumor qual
        // - excessive normal support
        if(belowHardMinQual(variant, genotypeIds, isSgl))
            return true;

        if(hasExcessiveReferenceSupport(variant, genotypeIds, isSgl))
            return true;

        return false;
    }

    public boolean belowHardMinQual(final VariantContext variant, final GenotypeIds genotypeIds, boolean isSgl)
    {
        Genotype tumorGenotype = variant.getGenotype(genotypeIds.TumorOrdinal);
        double qual = getGenotypeAttributeAsDouble(tumorGenotype, QUAL, 0);
        return qual < mFilterConstants.MinTumorQual;
    }

    private boolean hasExcessiveReferenceSupport(final VariantContext variant, final GenotypeIds genotypeIds, boolean isSgl)
    {
        if(!genotypeIds.hasReference() || mGermlineMode)
            return false;

        Genotype refGenotype = variant.getGenotype(genotypeIds.ReferenceOrdinal);
        Genotype tumorGenotype = variant.getGenotype(genotypeIds.TumorOrdinal);

        int refFrags = getGenotypeAttributeAsInt(refGenotype, TOTAL_FRAGS, 0);
        int tumorFrags = getGenotypeAttributeAsInt(tumorGenotype, TOTAL_FRAGS, 0);

        if(refFrags > mFilterConstants.HardMaxNormalRelativeSupport * tumorFrags)
            return true;

        if(refFrags > mFilterConstants.HardMaxNormalAbsoluteSupport && refFrags > mFilterConstants.SoftMaxNormalRelativeSupport * tumorFrags)
            return true;

        return false;
    }





}
