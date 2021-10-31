package com.hartwig.hmftools.svtools.germline;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
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
import static com.hartwig.hmftools.svtools.germline.VcfUtils.ASRP;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.HOMSEQ;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.IC;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.REFPAIR;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.SB;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.SR;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.getGenotypeAttributeAsInt;

import java.util.List;

import com.google.common.collect.Lists;

public class SoftFilters
{
    private final FilterConstants mFilterConstants;

    public SoftFilters(final FilterConstants filterConstants)
    {
        mFilterConstants = filterConstants;
    }

    public List<FilterType> applyFilters(final SvData sv)
    {
        List<FilterType> filters = Lists.newArrayList();

        if(normalCoverage(sv))
            filters.add(MIN_NORMAL_COVERAGE);

        if(normalRelativeSupport(sv))
            filters.add(MAX_NORMAL_RELATIVE_SUPPORT);

        if(allelicFrequency(sv))
            filters.add(MIN_TUMOR_AF);

        if(shortDelInsertArtifact(sv))
            filters.add(SHORT_DEL_INS_ARTIFACT);

        if(minQuality(sv))
            filters.add(MIN_QUAL);

        if(imprecise(sv))
            filters.add(IMPRECISE);

        if(polyGCInsert(sv))
            filters.add(MAX_POLY_G_LENGTH);

        if(polyATHomology(sv))
            filters.add(MAX_POLY_A_HOM_LENGTH);

        if(homologyLengthFilterShortInversion(sv))
            filters.add(MAX_HOM_LENGTH_SHORT_INV);

        if(inexactHomologyLengthShortDel(sv))
            filters.add(MAX_INEXACT_HOM_LENGTH_SHORT_DEL);

        if(shortSplitReadTumor(sv))
            filters.add(SHORT_SR_SUPPORT);

        if(shortSplitReadNormal(sv))
            filters.add(SHORT_SR_NORMAL);

        if(strandBias(sv))
            filters.add(SHORT_STRAND_BIAS);

        if(discordantPairSupport(sv) || breakendAssemblyReadPairs(sv))
            filters.add(DISCORDANT_PAIR_SUPPORT);

        if(minLength(sv))
            filters.add(MIN_LENGTH);

        return filters;
    }

    private boolean normalCoverage(final SvData sv)
    {
        if(!sv.hasReference())
            return false;

         return sv.referenceFragments() + sv.referenceReads() + sv.referencePairReads() < mFilterConstants.MinNormalCoverage;
    }

    private boolean normalRelativeSupport(final SvData sv)
    {
        return sv.referenceFragments() > mFilterConstants.SoftMaxNormalRelativeSupport * sv.tumorFragments();
    }

    private boolean allelicFrequency(final SvData sv)
    {
        int tumorFrags = sv.tumorFragments();
        int readPairSupport = (sv.isSgl() || !sv.isShortLocal()) ? sv.referencePairReads() : 0;
        int totalSupport = tumorFrags + sv.referenceReads() + readPairSupport;
        double alleleFrequency = totalSupport > 0 ? tumorFrags / totalSupport : 0;

        return alleleFrequency < mFilterConstants.MinTumorAf;
    }

    private boolean shortDelInsertArtifact(final SvData sv)
    {
        if(sv.type() != DEL || sv.length() >= SHORT_CALLING_SIZE)
            return false;

        return sv.altString().length() - 1 == sv.insertSequence().length();
    }

    private boolean minQuality(final SvData sv)
    {
        if(sv.isSgl())
            return sv.tumorQuality() < mFilterConstants.MinQualBreakend;
        else
            return sv.tumorQuality() < mFilterConstants.MinQualBreakpoint;
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
        String homologySequence = sv.contextStart().getAttributeAsString(HOMSEQ, "");
        return homologySequence.contains(POLY_A_HOMOLOGY) || homologySequence.contains(POLY_T_HOMOLOGY);
    }

    private boolean inexactHomologyLengthShortDel(final SvData sv)
    {

        /*
        fun inexactHomologyLengthShortDel(maxInexactHomLength: Int, minDelLength: Int = 100, maxDelLength: Int = 800): Boolean {
        return variantType is Deletion && variantType.length >= minDelLength && variantType.length <= maxDelLength && context.inexactHomologyLength() > maxInexactHomLength
        }

        fun inexactHomologyStart(): Int { return context.inexactHomologyStart(); }
        fun inexactHomologyEnd(): Int { return context.inexactHomologyEnd(); }

         */
        return false;
    }

    private boolean breakendAssemblyReadPairs(final SvData sv)
    {
        /*
                fun breakendAssemblyReadPairs(): Boolean {
        return isSingle && context.breakendAssemblyReadPairs() == 0 && !context.breakendAssemblyReadPairsIsInconsistent()
        }


         */

        return false;
    }

    private boolean imprecise(final SvData sv)
    {
        return sv.contextStart().getAttributeAsBoolean(VcfUtils.IMPRECISE, false);
    }

    private boolean homologyLengthFilterShortInversion(final SvData sv)
    {
        if(sv.type() != INV)
            return false;

        return sv.length() <= HOM_INV_LENGTH && sv.startHomology().length() > mFilterConstants.MaxHomLengthShortInv;
    }

    private boolean shortSplitReadTumor(final SvData sv)
    {
        int splitReads = getGenotypeAttributeAsInt(sv.tumorGenotype(), SR, 0);
        int indelCount = getGenotypeAttributeAsInt(sv.tumorGenotype(), IC, 0);

        return sv.isShortLocal() && (splitReads + indelCount == 0);
    }

    private boolean shortSplitReadNormal(final SvData sv)
    {
        if(!sv.hasReference())
            return false;

        int splitReads = getGenotypeAttributeAsInt(sv.refGenotype(), SR, 0);
        int indelCount = getGenotypeAttributeAsInt(sv.refGenotype(), IC, 0);

        return sv.isShortLocal() && (splitReads + indelCount > 0);
    }

    private boolean strandBias(final SvData sv)
    {
        if(sv.isShortLocal())
        {
            double strandBias = sv.contextStart().getAttributeAsDouble(SB, 0.5);
            return max(strandBias, 1 - strandBias) > mFilterConstants.MaxShortStrandBias;
        }

        return false;
    }

    private boolean discordantPairSupport(final SvData sv)
    {
        if(!sv.hasReference())
            return false;

        if(sv.isSgl() || sv.isShortLocal())
            return false;

        return sv.referencePairReads() == 0
                && VcfUtils.getGenotypeAttributeAsInt(sv.refGenotype(), ASRP, 0) == 0
                && VcfUtils.getGenotypeAttributeAsInt(sv.tumorGenotype(), REFPAIR, 0) == 0
                && VcfUtils.getGenotypeAttributeAsInt(sv.tumorGenotype(), ASRP, 0) == 0;
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
