package com.hartwig.hmftools.esvee.caller;

import static java.lang.Math.max;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.calcTrimmedBaseLength;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.INV_SHORT_LENGTH;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.INV_SHORT_MAX_HOMOLOGY;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.INV_SHORT_MIN_AF;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.MIN_AVG_FRAG_FACTOR;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.MIN_AVG_FRAG_STD_DEV_FACTOR;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.MIN_TRIMMED_ANCHOR_LENGTH;
import static com.hartwig.hmftools.esvee.common.FilterType.DUPLICATE;
import static com.hartwig.hmftools.esvee.common.FilterType.MIN_ANCHOR_LENGTH;
import static com.hartwig.hmftools.esvee.common.FilterType.MIN_LENGTH;
import static com.hartwig.hmftools.esvee.common.FilterType.MIN_QUALITY;
import static com.hartwig.hmftools.esvee.common.FilterType.MIN_AF;
import static com.hartwig.hmftools.esvee.common.FilterType.MIN_SUPPORT;
import static com.hartwig.hmftools.esvee.common.FilterType.SGL;
import static com.hartwig.hmftools.esvee.common.FilterType.SHORT_FRAG_LENGTH;
import static com.hartwig.hmftools.esvee.common.FilterType.SHORT_LOW_VAF_INV;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;
import com.hartwig.hmftools.esvee.common.FilterType;
import com.hartwig.hmftools.esvee.common.FragmentLengthBounds;

import htsjdk.variant.variantcontext.Genotype;

public class VariantFilters
{
    private final FilterConstants mFilterConstants;
    private final FragmentLengthBounds mFragmentLengthBounds;

    public VariantFilters(final FilterConstants filterConstants, final FragmentLengthBounds fragmentLengthBounds)
    {
        mFilterConstants = filterConstants;
        mFragmentLengthBounds = fragmentLengthBounds;
    }

    public void applyFilters(final Variant var)
    {
        // keep existing filters eg from assembly
        applyExistingFilters(var.breakendStart());

        if(!var.isSgl())
            applyExistingFilters(var.breakendEnd());

        if(mFilterConstants.FilterSGLs && var.isSgl())
            var.addFilter(SGL);

        if(belowMinAf(var))
            var.addFilter(MIN_AF);

        if(belowMinSupport(var))
            var.addFilter(MIN_SUPPORT);

        if(belowMinQuality(var))
            var.addFilter(MIN_QUALITY);

        if(belowMinLength(var))
            var.addFilter(MIN_LENGTH);

        if(belowMinAnchorLength(var))
            var.addFilter(MIN_ANCHOR_LENGTH);

        if(belowMinFragmentLength(var))
            var.addFilter(SHORT_FRAG_LENGTH);

        if(isShortLowVafInversion(var))
            var.addFilter(SHORT_LOW_VAF_INV);
    }

    private void applyExistingFilters(final Breakend breakend)
    {
        for(String filterStr : breakend.Context.getFilters())
        {
            try
            {
                FilterType filterType = FilterType.fromVcfTag(filterStr);

                if(filterType != null && filterType != DUPLICATE) // only required from where assembly marked duplicates
                    breakend.sv().addFilter(filterType);
            }
            catch(Exception e) {}
        }
    }

    private boolean belowMinSupport(final Variant var)
    {
        double supportThreshold;

        if(var.isHotspot())
        {
            supportThreshold = mFilterConstants.MinSupportHotspot;
        }
        else if(var.isSgl() && !var.isLineSite())
        {
            supportThreshold = mFilterConstants.MinSupportSgl;
        }
        else
        {
            supportThreshold = mFilterConstants.MinSupportJunction;
        }

        Breakend breakend = var.breakendStart();

        for(Genotype genotype : breakend.Context.getGenotypes())
        {
            int fragmentCount = breakend.fragmentCount(genotype);

            if(var.isSgl() && breakend.lineSiteBreakend() != null)
                fragmentCount += breakend.lineSiteBreakend().fragmentCount(genotype);

            if(fragmentCount >= supportThreshold)
                return false;
        }

        return true;
    }

    private boolean belowMinAf(final Variant var)
    {
        double afThreshold;

        if(var.isHotspot())
        {
            afThreshold = mFilterConstants.MinAfHotspot;
        }
        else if(var.isSgl() && !var.isLineSite())
        {
            afThreshold = mFilterConstants.MinAfSgl;
        }
        else
        {
            afThreshold = mFilterConstants.MinAfJunction;
        }

        return !anySamplesAboveAfThreshold(var, afThreshold);
    }

    private static boolean anySamplesAboveAfThreshold(final Variant var, final double afThreshold)
    {
        // both breakends for any sample must be above the threshold
        int genotypeCount = var.breakendStart().Context.getGenotypes().size();

        for(int g = 0; g < genotypeCount; ++g)
        {
            boolean aboveThreshold = true;

            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(var.breakends()[se] == null)
                    continue;

                Breakend breakend = var.breakends()[se];
                Genotype genotype = breakend.Context.getGenotypes().get(g);

                double af = breakend.calcAllelicFrequency(genotype);

                aboveThreshold &= af >= afThreshold;
            }

            if(aboveThreshold)
                return true;
        }

        return false;
    }

    private boolean belowMinQuality(final Variant var)
    {
        double qualThreshold = var.isHotspot() ? mFilterConstants.MinQualHotspot : mFilterConstants.MinQual;

        Breakend breakend = var.breakendStart();
        Breakend otherBreakend = var.breakendEnd();

        if(mFilterConstants.LowQualRegion.containsPosition(breakend.Chromosome, breakend.Position))
            qualThreshold *= 0.5;
        else if(otherBreakend != null && mFilterConstants.LowQualRegion.containsPosition(otherBreakend.Chromosome, otherBreakend.Position))
            qualThreshold *= 0.5;

        double variantQual = var.qual();

        if(breakend.lineSiteBreakend() != null)
            variantQual += breakend.lineSiteBreakend().sv().qual();

        return variantQual < qualThreshold;
    }

    private boolean belowMinLength(final Variant var)
    {
        if(var.type() == DEL || var.type() == DUP || var.type() == INS)
            return var.adjustedLength() < mFilterConstants.MinLength;

        return false;
    }

    private boolean belowMinAnchorLength(final Variant var)
    {
        if(var.isLineSite())
            return false;

        // skip for chained breakends
        if(var.inChainedAssembly())
            return false;

        if(belowMinAnchorLength(var.breakendStart()))
            return true;

        if(var.isSgl())
        {
            byte[] insertSequence = var.insertSequence().getBytes();
            List<RepeatInfo> repeats = RepeatInfo.findRepeats(insertSequence);
            int trimmedSequenceLength = calcTrimmedBaseLength(0, insertSequence.length - 1, repeats);

            return trimmedSequenceLength < MIN_TRIMMED_ANCHOR_LENGTH;
        }
        else
        {
            return belowMinAnchorLength(var.breakendEnd());
        }
    }

    private boolean belowMinAnchorLength(final Breakend breakend)
    {
        return breakend.anchorLength() < MIN_TRIMMED_ANCHOR_LENGTH;
    }

    private boolean belowMinFragmentLength(final Variant var)
    {
        if(var.isSgl())
            return false;

        if(var.isLineSite())
            return false;

        if(var.inChainedAssembly())
            return false;

        int svAvgLength = var.averageFragmentLength();

        if(svAvgLength == 0) // for now treat is this as a pass
            return false;

        int medianLength = mFragmentLengthBounds.Median;
        double stdDeviation = mFragmentLengthBounds.StdDeviation;

        int totalSplitFrags = var.splitFragmentCount();
        double lowerStdDevThreshold = MIN_AVG_FRAG_STD_DEV_FACTOR * stdDeviation;
        double lowerLengthLimit = medianLength - max(MIN_AVG_FRAG_FACTOR * stdDeviation / sqrt(totalSplitFrags), lowerStdDevThreshold);

        return svAvgLength < lowerLengthLimit;
    }

    private boolean isShortLowVafInversion(final Variant var)
    {
        // FILTER if [TYPE=INV, AF<0.05, LEN<2500
        if(var.type() != INV)
            return false;

        if(var.adjustedLength() > INV_SHORT_LENGTH)
            return false;

        if(anySamplesAboveAfThreshold(var, INV_SHORT_MIN_AF))
            return false;

        int inexactHomology = var.breakendStart().InexactHomology.length();
        return inexactHomology > INV_SHORT_MAX_HOMOLOGY;
    }

    public static void logFilterTypeCounts(final List<Variant> variantList)
    {
        Map<FilterType,Integer> filterCounts = Maps.newHashMap();

        for(Variant var : variantList)
        {
            for(FilterType filterType : var.filters())
            {
                int count = filterCounts.getOrDefault(filterType, 0);
                filterCounts.put(filterType, count + 1);
            }
        }

        for(Map.Entry<FilterType,Integer> entry : filterCounts.entrySet())
        {
            SV_LOGGER.debug("variant filter {}: count({})", entry.getKey(), entry.getValue());
        }
    }
}
