package com.hartwig.hmftools.svtools.germline;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.svtools.germline.FilterConstants.HOM_INV_LENGTH;
import static com.hartwig.hmftools.svtools.germline.FilterConstants.POLY_A_HOMOLOGY;
import static com.hartwig.hmftools.svtools.germline.FilterConstants.POLY_C_INSERT;
import static com.hartwig.hmftools.svtools.germline.FilterConstants.POLY_G_INSERT;
import static com.hartwig.hmftools.svtools.germline.FilterConstants.POLY_T_HOMOLOGY;
import static com.hartwig.hmftools.svtools.germline.FilterConstants.SHORT_CALLING_SIZE;
import static com.hartwig.hmftools.svtools.germline.FilterType.DISCORDANT_PAIR_SUPPORT;
import static com.hartwig.hmftools.svtools.germline.FilterType.IMPRECISE;
import static com.hartwig.hmftools.svtools.germline.FilterType.MAX_HOM_LENGTH_SHORT_INV;
import static com.hartwig.hmftools.svtools.germline.FilterType.MAX_INEXACT_HOM_LENGTH_SHORT_DEL;
import static com.hartwig.hmftools.svtools.germline.FilterType.MAX_NORMAL_RELATIVE_SUPPORT;
import static com.hartwig.hmftools.svtools.germline.FilterType.MAX_POLY_A_HOM_LENGTH;
import static com.hartwig.hmftools.svtools.germline.FilterType.MAX_POLY_G_LENGTH;
import static com.hartwig.hmftools.svtools.germline.FilterType.MIN_LENGTH;
import static com.hartwig.hmftools.svtools.germline.FilterType.MIN_NORMAL_COVERAGE;
import static com.hartwig.hmftools.svtools.germline.FilterType.MIN_QUAL;
import static com.hartwig.hmftools.svtools.germline.FilterType.MIN_TUMOR_AF;
import static com.hartwig.hmftools.svtools.germline.FilterType.SHORT_DEL_INS_ARTIFACT;
import static com.hartwig.hmftools.svtools.germline.FilterType.SHORT_SR_NORMAL;
import static com.hartwig.hmftools.svtools.germline.FilterType.SHORT_SR_SUPPORT;
import static com.hartwig.hmftools.svtools.germline.FilterType.SHORT_STRAND_BIAS;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.VT_ASRP;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.VT_HOMSEQ;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.VT_IC;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.VT_REFPAIR;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.VT_SB;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.VT_SR;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.getGenotypeAttributeAsInt;

public class SoftFilters
{
    private final FilterConstants mFilterConstants;

    public SoftFilters(final FilterConstants filterConstants)
    {
        mFilterConstants = filterConstants;
    }

    public void applyFilters(final SvData sv)
    {
        // breakend filters
        for(int se = SE_START; se <= SE_END; ++se)
        {
            Breakend breakend = sv.breakends()[se];

            if(normalCoverage(breakend))
                breakend.addFilter(MIN_NORMAL_COVERAGE);

            if(normalRelativeSupport(breakend))
                breakend.addFilter(MAX_NORMAL_RELATIVE_SUPPORT);

            if(allelicFrequency(sv, breakend))
                breakend.addFilter(MIN_TUMOR_AF);

            if(minQuality(sv, breakend))
                breakend.addFilter(MIN_QUAL);

            if(shortSplitReadTumor(sv, breakend))
                breakend.addFilter(SHORT_SR_SUPPORT);

            if(shortSplitReadNormal(sv, breakend))
                breakend.addFilter(SHORT_SR_NORMAL);

            if(discordantPairSupport(sv, breakend) || breakendAssemblyReadPairs(breakend))
                breakend.addFilter(DISCORDANT_PAIR_SUPPORT);

            if(shortDelInsertArtifact(sv, breakend))
                sv.addFilter(SHORT_DEL_INS_ARTIFACT);
        }

        // SV filters
        if(imprecise(sv))
            sv.addFilter(IMPRECISE);

        if(polyGCInsert(sv))
            sv.addFilter(MAX_POLY_G_LENGTH);

        if(polyATHomology(sv))
            sv.addFilter(MAX_POLY_A_HOM_LENGTH);

        if(homologyLengthFilterShortInversion(sv))
            sv.addFilter(MAX_HOM_LENGTH_SHORT_INV);

        if(inexactHomologyLengthShortDel(sv))
            sv.addFilter(MAX_INEXACT_HOM_LENGTH_SHORT_DEL);

        if(strandBias(sv))
            sv.addFilter(SHORT_STRAND_BIAS);

        if(minLength(sv))
            sv.addFilter(MIN_LENGTH);
    }

    private boolean normalCoverage(final Breakend breakend)
    {
        if(breakend.RefGenotype == null)
            return false;

         return breakend.ReferenceFragments + breakend.ReferenceReads + breakend.ReferencePairReads < mFilterConstants.MinNormalCoverage;
    }

    private boolean normalRelativeSupport(final Breakend breakend)
    {
        if(breakend.RefGenotype == null)
            return false;

        return breakend.ReferenceFragments > mFilterConstants.SoftMaxNormalRelativeSupport * breakend.TumorFragments;
    }

    private boolean allelicFrequency(final SvData sv, final Breakend breakend)
    {
        int tumorFrags = breakend.TumorFragments;
        int readPairSupport = (sv.isSgl() || !sv.isShortLocal()) ? breakend.ReferencePairReads : 0;
        int totalSupport = tumorFrags + breakend.ReferenceReads + readPairSupport;
        double alleleFrequency = totalSupport > 0 ? tumorFrags / totalSupport : 0;

        return alleleFrequency < mFilterConstants.MinTumorAf;
    }

    private boolean shortDelInsertArtifact(final SvData sv, final Breakend breakend)
    {
        if(sv.type() != DEL || sv.length() >= SHORT_CALLING_SIZE)
            return false;

        return breakend.Alt.length() - 1 == breakend.insertSequenceLength();
    }

    private boolean minQuality(final SvData sv, final Breakend breakend)
    {
        if(sv.isSgl())
            return breakend.Qual < mFilterConstants.MinQualBreakend;
        else
            return breakend.Qual < mFilterConstants.MinQualBreakpoint;
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
            if(mFilterConstants.PolyGcRegion.containsPosition(sv.chromosomeEnd(), sv.posEnd()));
            return true;
        }

        return false;
    }

    private boolean polyATHomology(final SvData sv)
    {
        String homologySequence = sv.contextStart().getAttributeAsString(VT_HOMSEQ, "");
        return homologySequence.contains(POLY_A_HOMOLOGY) || homologySequence.contains(POLY_T_HOMOLOGY);
    }

    private boolean inexactHomologyLengthShortDel(final SvData sv)
    {
        // TODO

        /*
        fun inexactHomologyLengthShortDel(maxInexactHomLength: Int, minDelLength: Int = 100, maxDelLength: Int = 800): Boolean {
        return variantType is Deletion && variantType.length >= minDelLength && variantType.length <= maxDelLength && context.inexactHomologyLength() > maxInexactHomLength
        }

        fun inexactHomologyStart(): Int { return context.inexactHomologyStart(); }
        fun inexactHomologyEnd(): Int { return context.inexactHomologyEnd(); }

         */
        return false;
    }

    private boolean breakendAssemblyReadPairs(final Breakend breakend)
    {
        // TODO
        /*
                fun breakendAssemblyReadPairs(): Boolean {
        return isSingle && context.breakendAssemblyReadPairs() == 0 && !context.breakendAssemblyReadPairsIsInconsistent()
        }


         */

        return false;
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

    private boolean strandBias(final SvData sv)
    {
        if(sv.isShortLocal())
        {
            double strandBias = sv.contextStart().getAttributeAsDouble(VT_SB, 0.5);
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

        return breakend.ReferencePairReads == 0
                && VcfUtils.getGenotypeAttributeAsInt(breakend.RefGenotype, VT_ASRP, 0) == 0
                && VcfUtils.getGenotypeAttributeAsInt(breakend.TumorGenotype, VT_REFPAIR, 0) == 0
                && VcfUtils.getGenotypeAttributeAsInt(breakend.TumorGenotype, VT_ASRP, 0) == 0;
    }

    private boolean minLength(final SvData sv)
    {
        if(sv.type() == DEL)
            return sv.length() + sv.insertSequence().length() - 1 < mFilterConstants.MinLength;
        else if(sv.type() == DUP)
            return sv.length() + sv.insertSequence().length() + 1 < mFilterConstants.MinLength;
        else if(sv.type() == INS)
            return sv.length() + sv.insertSequence().length() < mFilterConstants.MinLength;
        else
            return false;
    }
}
