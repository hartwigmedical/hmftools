package com.hartwig.hmftools.esvee.caller;

import static java.lang.Math.max;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.STRAND_BIAS;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FilterType.MIN_LENGTH;
import static com.hartwig.hmftools.esvee.common.FilterType.MIN_QUALITY;
import static com.hartwig.hmftools.esvee.common.FilterType.MIN_AF;
import static com.hartwig.hmftools.esvee.common.FilterType.MIN_SUPPORT;
import static com.hartwig.hmftools.esvee.common.FilterType.SGL;
import static com.hartwig.hmftools.esvee.common.FilterType.SHORT_FRAG_LENGTH;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.common.FilterType;
import com.hartwig.hmftools.esvee.common.FragmentLengthBounds;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

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
        if(mFilterConstants.FilterSGLs && var.isSgl())
            var.addFilter(SGL);

        if(belowMinAf(var))
            var.addFilter(MIN_AF);

        if(belowMinSupport(var))
            var.addFilter(MIN_SUPPORT);

        if(belowMinQuality(var))
            var.addFilter(MIN_QUALITY);

        if(hasStrandBias(var))
            var.addFilter(FilterType.STRAND_BIAS);

        if(belowMinLength(var))
            var.addFilter(MIN_LENGTH);

        if(belowMinFragmentLength(var))
            var.addFilter(SHORT_FRAG_LENGTH);
    }

    private boolean belowMinSupport(final Variant var)
    {
        double supportThreshold = var.isHotspot() ? mFilterConstants.MinSupportHotspot :
                (var.isSgl() ? mFilterConstants.MinSupportSgl : mFilterConstants.MinSupportJunction);

        Breakend breakend = var.breakendStart();

        for(Genotype genotype : breakend.Context.getGenotypes())
        {
            int fragmentCount = breakend.fragmentCount(genotype);

            if(fragmentCount >= supportThreshold)
                return false;
        }

        return true;
    }

    private boolean belowMinAf(final Variant var)
    {
        double afThreshold = var.isHotspot() ? mFilterConstants.MinAfHotspot :
                (var.isSgl() ? mFilterConstants.MinAfSgl : mFilterConstants.MinAfJunction);

        Breakend breakend = var.breakendStart();

        for(Genotype genotype : breakend.Context.getGenotypes())
        {
            double af = breakend.calcAllelicFrequency(genotype);

            if(af >= afThreshold)
                return false;
        }

        return true;
    }

    private boolean belowMinQuality(final Variant var)
    {
        double qualThreshold = mFilterConstants.MinQual;

        Breakend breakend = var.breakendStart();
        Breakend otherBreakend = var.breakendEnd();

        if(mFilterConstants.LowQualRegion.containsPosition(breakend.Chromosome, breakend.Position))
            qualThreshold *= 0.5;
        else if(otherBreakend != null && mFilterConstants.LowQualRegion.containsPosition(otherBreakend.Chromosome, otherBreakend.Position))
            qualThreshold *= 0.5;

        return var.qual() < qualThreshold;
    }

    private boolean belowMinLength(final Variant var)
    {
        if(var.type() == DEL)
            return var.length() + var.insertSequence().length() - 1 < mFilterConstants.MinLength;
        else if(var.type() == DUP)
            return var.length() + var.insertSequence().length() < mFilterConstants.MinLength;
        else if(var.type() == INS)
            return var.length() + var.insertSequence().length() + 1 < mFilterConstants.MinLength;
        else
            return false;
    }

    private boolean belowMinFragmentLength(final Variant var)
    {
        int medianLength = mFragmentLengthBounds.Median;
        double stdDeviation = mFragmentLengthBounds.StdDeviation;

        int totalSplitFrags = var.splitFragmentCount();
        double lowerLengthLimit = medianLength - (mFilterConstants.MinAvgFragFactor * stdDeviation / sqrt(totalSplitFrags));

        int svAvgLength = var.averageFragmentLength();

        return svAvgLength < lowerLengthLimit;
    }

    private boolean hasStrandBias(final Variant var)
    {
        /*
        private boolean singleStrandBias(final Breakend breakend)
        {
            if(!breakend.isSgl() || breakend.IsLineInsertion)
                return false;

            if(mFilterConstants.LowQualRegion.containsPosition(breakend.Chromosome, breakend.Position))
                return false;

            double strandBias = calcStrandBias(breakend.Context);
            return strandBias < SGL_MIN_STRAND_BIAS || strandBias > SGL_MAX_STRAND_BIAS;
        }

        if(sv.isShortLocal())
        {
            double strandBias = breakend.Context.getAttributeAsDouble(STRAND_BIAS, 0.5);
            return max(strandBias, 1 - strandBias) > MAX_STRAND_BIAS;
        }
        */

        return false;
    }

    private static double calcStrandBias(final VariantContext variantContext)
    {
        double strandBias = variantContext.getAttributeAsDouble(STRAND_BIAS, 0.5);
        return max(strandBias, 1 - strandBias);
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
