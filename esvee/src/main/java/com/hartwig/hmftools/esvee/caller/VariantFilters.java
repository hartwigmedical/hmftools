package com.hartwig.hmftools.esvee.caller;

import static java.lang.Math.max;

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

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.common.FilterType;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class VariantFilters
{
    private final FilterConstants mFilterConstants;

    public VariantFilters(final FilterConstants filterConstants)
    {
        mFilterConstants = filterConstants;
    }

    public void applyFilters(final SvData sv)
    {
        if(mFilterConstants.FilterSGLs && sv.isSgl())
            sv.addFilter(SGL);

        if(belowMinAf(sv))
            sv.addFilter(MIN_AF);

        if(belowMinSupport(sv))
            sv.addFilter(MIN_SUPPORT);

        if(belowMinQuality(sv))
            sv.addFilter(MIN_QUALITY);

        if(hasStrandBias(sv))
            sv.addFilter(FilterType.STRAND_BIAS);

        if(belowMinLength(sv))
            sv.addFilter(MIN_LENGTH);
    }

    private boolean belowMinSupport(final SvData sv)
    {
        double supportThreshold = sv.isHotspot() ? mFilterConstants.MinSupportHotspot :
                (sv.isSgl() ? mFilterConstants.MinSupportSgl : mFilterConstants.MinSupportJunction);

        Breakend breakend = sv.breakendStart();

        for(Genotype genotype : breakend.Context.getGenotypes())
        {
            int fragmentCount = breakend.fragmentCount(genotype);

            if(fragmentCount >= supportThreshold)
                return false;
        }

        return true;
    }

    private boolean belowMinAf(final SvData sv)
    {
        double afThreshold = sv.isHotspot() ? mFilterConstants.MinAfHotspot :
                (sv.isSgl() ? mFilterConstants.MinAfSgl : mFilterConstants.MinAfJunction);

        Breakend breakend = sv.breakendStart();

        for(Genotype genotype : breakend.Context.getGenotypes())
        {
            double af = breakend.calcAllelicFrequency(genotype);

            if(af >= afThreshold)
                return false;
        }

        return true;
    }

    private boolean belowMinQuality(final SvData sv)
    {
        double qualThreshold = mFilterConstants.MinQual;

        Breakend breakend = sv.breakendStart();
        Breakend otherBreakend = sv.breakendEnd();

        if(mFilterConstants.LowQualRegion.containsPosition(breakend.Chromosome, breakend.Position))
            qualThreshold *= 0.5;
        else if(otherBreakend != null && mFilterConstants.LowQualRegion.containsPosition(otherBreakend.Chromosome, otherBreakend.Position))
            qualThreshold *= 0.5;

        return sv.qual() < qualThreshold;
    }

    private boolean belowMinLength(final SvData sv)
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

    private boolean hasStrandBias(final SvData sv)
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

    public static void logFilterTypeCounts(final List<SvData> svDataList)
    {
        Map<FilterType,Integer> filterCounts = Maps.newHashMap();

        for(SvData var : svDataList)
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
