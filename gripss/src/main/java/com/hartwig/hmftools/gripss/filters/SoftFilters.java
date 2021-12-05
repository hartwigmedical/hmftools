package com.hartwig.hmftools.gripss.filters;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_REF;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_RP;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.HOM_INV_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.INEXACT_HOM_LENGTH_SHORT_DEL_MAX_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.INEXACT_HOM_LENGTH_SHORT_DEL_MIN_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.POLY_A_HOMOLOGY;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.POLY_C_INSERT;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.POLY_G_INSERT;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.POLY_T_HOMOLOGY;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.SGL_INS_SEQ_MIN_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.SGL_MAX_STRAND_BIAS;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.SGL_MIN_STRAND_BIAS;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.SHORT_CALLING_SIZE;
import static com.hartwig.hmftools.gripss.filters.FilterType.DISCORDANT_PAIR_SUPPORT;
import static com.hartwig.hmftools.gripss.filters.FilterType.IMPRECISE;
import static com.hartwig.hmftools.gripss.filters.FilterType.MAX_HOM_LENGTH_SHORT_INV;
import static com.hartwig.hmftools.gripss.filters.FilterType.MAX_INEXACT_HOM_LENGTH_SHORT_DEL;
import static com.hartwig.hmftools.gripss.filters.FilterType.MAX_NORMAL_RELATIVE_SUPPORT;
import static com.hartwig.hmftools.gripss.filters.FilterType.MAX_POLY_A_HOM_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterType.MAX_POLY_G_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterType.MIN_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterType.MIN_NORMAL_COVERAGE;
import static com.hartwig.hmftools.gripss.filters.FilterType.MIN_QUAL;
import static com.hartwig.hmftools.gripss.filters.FilterType.MIN_TUMOR_AF;
import static com.hartwig.hmftools.gripss.filters.FilterType.SGL_INSERT_SEQ_MIN_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterType.SGL_STRAND_BIAS;
import static com.hartwig.hmftools.gripss.filters.FilterType.SHORT_DEL_INS_ARTIFACT;
import static com.hartwig.hmftools.gripss.filters.FilterType.SHORT_SR_NORMAL;
import static com.hartwig.hmftools.gripss.filters.FilterType.SHORT_SR_SUPPORT;
import static com.hartwig.hmftools.gripss.filters.FilterType.SHORT_STRAND_BIAS;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_ASRP;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_HOMSEQ;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_IC;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_REFPAIR;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_SB;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_SR;
import static com.hartwig.hmftools.gripss.common.VcfUtils.getGenotypeAttributeAsInt;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.gripss.FilterCache;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.common.VcfUtils;

import htsjdk.variant.variantcontext.VariantContext;

public class SoftFilters
{
    private final FilterConstants mFilterConstants;

    public SoftFilters(final FilterConstants filterConstants)
    {
        mFilterConstants = filterConstants;
    }

    public void applyFilters(final SvData sv, final FilterCache filterCache)
    {
        if(filterCache.isHotspot(sv))
            return;

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
                filters.add(MIN_QUAL);

            if(shortSplitReadTumor(sv, breakend))
                filters.add(SHORT_SR_SUPPORT);

            if(shortSplitReadNormal(sv, breakend))
                filters.add(SHORT_SR_NORMAL);

            if(discordantPairSupport(sv, breakend))
                filters.add(DISCORDANT_PAIR_SUPPORT);

            if(singleStrandBias(breakend))
                filters.add(SGL_STRAND_BIAS);

            if(singleInsertSequenceMinLength(breakend))
                filters.add(SGL_INSERT_SEQ_MIN_LENGTH);

            if(shortDelInsertArtifact(sv, breakend))
                filters.add(SHORT_DEL_INS_ARTIFACT);

            if(inexactHomologyLengthShortDel(sv, breakend))
                filters.add(MAX_INEXACT_HOM_LENGTH_SHORT_DEL);

            if(strandBias(sv, breakend))
                filters.add(SHORT_STRAND_BIAS);

            if((se == SE_END && beStartFilters.contains(IMPRECISE)) || imprecise(sv))
                filters.add(IMPRECISE);

            if((se == SE_END && beStartFilters.contains(MAX_POLY_G_LENGTH)) || polyGCInsert(sv))
                filters.add(MAX_POLY_G_LENGTH);

            if((se == SE_END && beStartFilters.contains(MAX_POLY_A_HOM_LENGTH)) || polyATHomology(sv))
                filters.add(MAX_POLY_A_HOM_LENGTH);

            if((se == SE_END && beStartFilters.contains(MAX_HOM_LENGTH_SHORT_INV)) || homologyLengthFilterShortInversion(sv))
                filters.add(MAX_HOM_LENGTH_SHORT_INV);

            if((se == SE_END && beStartFilters.contains(MIN_LENGTH)) || minLength(sv))
                filters.add(MIN_LENGTH);

            if(!filters.isEmpty())
                filterCache.addBreakendFilters(breakend, filters);
        }
    }

    private boolean normalCoverage(final Breakend breakend)
    {
        if(breakend.RefGenotype == null)
            return false;

        int refSupportReads = getGenotypeAttributeAsInt(breakend.RefGenotype, VT_REF, 0);
        int refSupportReadPairs = getGenotypeAttributeAsInt(breakend.RefGenotype, VT_REFPAIR, 0);

        return breakend.ReferenceFragments + refSupportReads + refSupportReadPairs < mFilterConstants.MinNormalCoverage;
    }

    private boolean normalRelativeSupport(final Breakend breakend)
    {
        if(breakend.RefGenotype == null)
            return false;

        return breakend.ReferenceFragments > mFilterConstants.SoftMaxNormalRelativeSupport * breakend.TumorFragments;
    }

    private boolean allelicFrequency(final SvData sv, final Breakend breakend)
    {
        double afThreshold = sv.isSgl() ? mFilterConstants.MinTumorAfBreakend : mFilterConstants.MinTumorAfBreakpoint;
        return breakend.AllelicFrequency < afThreshold;
    }

    private boolean shortDelInsertArtifact(final SvData sv, final Breakend breakend)
    {
        if(sv.type() != DEL)
            return false;

        int length = sv.length(); // lengths of 1 were treated as INS in gripsKT even without an insert sequence
        return length < SHORT_CALLING_SIZE && length > 1 && (length - 1 == breakend.insertSequenceLength());
    }

    private boolean minQuality(final SvData sv, final Breakend breakend)
    {
        double qualThreshold = sv.isSgl() ? mFilterConstants.MinQualBreakend : mFilterConstants.MinQualBreakpoint;

        if(mFilterConstants.LowQualRegion.containsPosition(breakend.Chromosome, breakend.Position))
            qualThreshold *= 0.5;

        return breakend.Qual < qualThreshold;
    }

    private boolean polyGCInsert(final SvData sv)
    {
        if(mFilterConstants.PolyGcRegion.containsPosition(sv.chromosomeStart(), sv.posStart()))
            return true;

        if(sv.isSgl())
        {
            if(sv.insertSequence().contains(POLY_G_INSERT) || sv.insertSequence().contains(POLY_C_INSERT))
                return true;
        }
        else
        {
            if(mFilterConstants.PolyGcRegion.containsPosition(sv.chromosomeStart(), sv.posStart())
            || mFilterConstants.PolyGcRegion.containsPosition(sv.chromosomeEnd(), sv.posEnd()))
            {
                return true;
            }
        }

        return false;
    }

    private boolean polyATHomology(final SvData sv)
    {
        String homologySequence = sv.contextStart().getAttributeAsString(VT_HOMSEQ, "");
        return homologySequence.contains(POLY_A_HOMOLOGY) || homologySequence.contains(POLY_T_HOMOLOGY);
    }

    private boolean inexactHomologyLengthShortDel(final SvData sv, final Breakend breakend)
    {
        if(sv.type() != DEL)
            return false;

        if(sv.length() < INEXACT_HOM_LENGTH_SHORT_DEL_MIN_LENGTH || sv.length() > INEXACT_HOM_LENGTH_SHORT_DEL_MAX_LENGTH)
            return false;

        return breakend.inexactHomologyLength() > mFilterConstants.MaxInexactHomLengthShortDel;
    }

    private static double calcStrandBias(final VariantContext variantContext)
    {
        double strandBias = variantContext.getAttributeAsDouble(VT_SB, 0.5);
        return max(strandBias, 1 - strandBias);
    }

    private boolean singleStrandBias(final Breakend breakend)
    {
        if(!breakend.isSgl() || breakend.IsLineInsertion)
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

    private boolean imprecise(final SvData sv)
    {
        return sv.imprecise();
    }

    private boolean homologyLengthFilterShortInversion(final SvData sv)
    {
        if(sv.type() != INV)
            return false;

        return sv.length() <= HOM_INV_LENGTH && sv.startHomology().length() > mFilterConstants.MaxHomLengthShortInv;
    }

    private boolean shortSplitReadTumor(final SvData sv, final Breakend breakend)
    {
        int splitReads = getGenotypeAttributeAsInt(breakend.TumorGenotype, VT_SR, 0);
        int indelCount = getGenotypeAttributeAsInt(breakend.TumorGenotype, VT_IC, 0);

        return sv.isShortLocal() && (splitReads + indelCount == 0);
    }

    private boolean shortSplitReadNormal(final SvData sv, final Breakend breakend)
    {
        if(breakend.RefGenotype == null)
            return false;

        int splitReads = getGenotypeAttributeAsInt(breakend.RefGenotype, VT_SR, 0);
        int indelCount = getGenotypeAttributeAsInt(breakend.RefGenotype, VT_IC, 0);

        return sv.isShortLocal() && (splitReads + indelCount > 0);
    }

    private boolean strandBias(final SvData sv, final Breakend breakend)
    {
        if(sv.isShortLocal())
        {
            double strandBias = breakend.Context.getAttributeAsDouble(VT_SB, 0.5);
            return max(strandBias, 1 - strandBias) > mFilterConstants.MaxShortStrandBias;
        }

        return false;
    }

    private boolean discordantPairSupport(final SvData sv, final Breakend breakend)
    {
        if(!sv.hasReference())
            return false;

        if(sv.isSgl() || sv.isShortLocal())
            return false;

        return VcfUtils.getGenotypeAttributeAsInt(breakend.RefGenotype, VT_RP, 0) == 0
                && VcfUtils.getGenotypeAttributeAsInt(breakend.RefGenotype, VT_ASRP, 0) == 0
                && VcfUtils.getGenotypeAttributeAsInt(breakend.TumorGenotype, VT_RP, 0) == 0
                && VcfUtils.getGenotypeAttributeAsInt(breakend.TumorGenotype, VT_ASRP, 0) == 0;
    }

    private boolean minLength(final SvData sv)
    {
        if(sv.type() == DEL)
            return sv.length() + sv.insertSequence().length() - 1 < mFilterConstants.MinLength;
        else if(sv.type() == DUP)
            return sv.length() + sv.insertSequence().length()  < mFilterConstants.MinLength;
        else if(sv.type() == INS)
            return sv.length() + sv.insertSequence().length() + 1 < mFilterConstants.MinLength;
        else
            return false;
    }
}
