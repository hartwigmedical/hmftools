package com.hartwig.hmftools.gripss.filters;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.ExcludedRegions.POLY_C_INSERT;
import static com.hartwig.hmftools.common.region.ExcludedRegions.POLY_G_INSERT;
import static com.hartwig.hmftools.common.sv.LineElements.POLY_A_HOMOLOGY;
import static com.hartwig.hmftools.common.sv.LineElements.POLY_T_HOMOLOGY;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOMSEQ;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.STRAND_BIAS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_ASRP;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_ASSR;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.INDEL_COUNT;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.READ_PAIRS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.SPLIT_READS;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.HOM_INV_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.SGL_INS_SEQ_MIN_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.SGL_MAX_STRAND_BIAS;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.SGL_MIN_STRAND_BIAS;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.SHORT_CALLING_SIZE;
import static com.hartwig.hmftools.gripss.filters.FilterType.DISCORDANT_PAIR_SUPPORT;
import static com.hartwig.hmftools.gripss.filters.FilterType.IMPRECISE;
import static com.hartwig.hmftools.gripss.filters.FilterType.MAX_HOM_LENGTH_SHORT_INV;
import static com.hartwig.hmftools.gripss.filters.FilterType.MAX_NORMAL_RELATIVE_SUPPORT;
import static com.hartwig.hmftools.gripss.filters.FilterType.MAX_POLY_A_HOM_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterType.MAX_POLY_G_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterType.MIN_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterType.MIN_NORMAL_COVERAGE;
import static com.hartwig.hmftools.gripss.filters.FilterType.MIN_QUAL;
import static com.hartwig.hmftools.gripss.filters.FilterType.MIN_TUMOR_AF;
import static com.hartwig.hmftools.gripss.filters.FilterType.MODIFIED_AF;
import static com.hartwig.hmftools.gripss.filters.FilterType.QUAL_PER_AD;
import static com.hartwig.hmftools.gripss.filters.FilterType.SGL_INSERT_SEQ_MIN_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterType.SGL_STRAND_BIAS;
import static com.hartwig.hmftools.gripss.filters.FilterType.SHORT_DEL_INS_ARTIFACT;
import static com.hartwig.hmftools.gripss.filters.FilterType.SHORT_SR_NORMAL;
import static com.hartwig.hmftools.gripss.filters.FilterType.SHORT_SR_SUPPORT;
import static com.hartwig.hmftools.gripss.filters.FilterType.SHORT_STRAND_BIAS;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.gripss.FilterCache;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.SvData;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class SoftFilters
{
    private final FilterConstants mFilterConstants;
    private final boolean mGermlineMode;

    public SoftFilters(final FilterConstants filterConstants, final boolean germlineMode)
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

            if(strandBias(sv, breakend))
                filters.add(SHORT_STRAND_BIAS);

            // the following filters are replicated in the end breakend if met in the start
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

            if((se == SE_END && beStartFilters.contains(MODIFIED_AF)) || modifiedAF(sv, false))
                filters.add(MODIFIED_AF);

            if(qualPerAD(breakend))
                filters.add(QUAL_PER_AD);

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

    private boolean qualPerAD(final Breakend breakend)
    {
        if(mFilterConstants.QualPerAD == 0)
            return false;

        int indelCount = getGenotypeAttributeAsInt(breakend.TumorGenotype, INDEL_COUNT, 0);
        int ad = indelCount + breakend.TumorFragments;

        if(ad == 0)
            return false;

        double qualPerAD = breakend.Qual / ad;
        return qualPerAD < mFilterConstants.QualPerAD;
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
        int indelCount = getGenotypeAttributeAsInt(breakend.TumorGenotype, INDEL_COUNT, 0);
        int ad = indelCount + breakend.TumorFragments;
        double totalSupport = ad + refSupport + readPairSupport;

        return totalSupport > 0 ? ad / totalSupport : 0;
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
        if(mFilterConstants.matchesPolyGRegion(sv.chromosomeStart(), sv.posStart()))
            return true;

        if(sv.isSgl())
        {
            if(sv.insertSequence().contains(POLY_G_INSERT) || sv.insertSequence().contains(POLY_C_INSERT))
                return true;
        }
        else
        {
            if(mFilterConstants.matchesPolyGRegion(sv.chromosomeEnd(), sv.posEnd()))
                return true;
        }

        return false;
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
        return sv.isShortLocal() && getSplitReadCount(breakend, breakend.TumorGenotype) == 0;
    }

    private boolean shortSplitReadNormal(final SvData sv, final Breakend breakend)
    {
        if(breakend.RefGenotype == null || mGermlineMode)
            return false;

        return sv.isShortLocal() && getSplitReadCount(breakend, breakend.RefGenotype) > 0;
    }

    private static int getSplitReadCount(final Breakend breakend, final Genotype genotype)
    {
        int splitReads = getGenotypeAttributeAsInt(genotype, SPLIT_READS, 0);
        int assemblySplitReads = getGenotypeAttributeAsInt(genotype, GRIDSS_ASSR, 0);
        int indelCount = getGenotypeAttributeAsInt(genotype, INDEL_COUNT, 0);
        return splitReads + assemblySplitReads + indelCount;
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

    private boolean discordantPairSupport(final SvData sv, final Breakend breakend)
    {
        if(!sv.hasReference() || mGermlineMode)
            return false;

        if(sv.type() != INV || sv.length() > HOM_INV_LENGTH)
            return false;

        return getGenotypeAttributeAsInt(breakend.RefGenotype, READ_PAIRS, 0) == 0
                && getGenotypeAttributeAsInt(breakend.RefGenotype, GRIDSS_ASRP, 0) == 0
                && getGenotypeAttributeAsInt(breakend.TumorGenotype, READ_PAIRS, 0) == 0
                && getGenotypeAttributeAsInt(breakend.TumorGenotype, GRIDSS_ASRP, 0) == 0;
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
}
