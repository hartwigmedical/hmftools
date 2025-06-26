package com.hartwig.hmftools.esvee.caller;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsDouble;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.calcTrimmedBaseLength;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.DEL_ARTEFACT_LENGTH_FACTOR;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.DEL_ARTEFACT_MAX_HOMOLOGY;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.DEL_ARTEFACT_MIN_AF;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.DEL_ARTEFACT_SHORT_LENGTH;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.INV_ADJACENT_EXCLUDED_FILTERS;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.INV_ADJACENT_LENGTH;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.INV_ADJACENT_MIN_UPS;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.INV_SHORT_FRAGMENT_AF_RATIO;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.INV_SHORT_FRAGMENT_LENGTH;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.INV_SHORT_FRAGMENT_MIN_AF;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.INV_SHORT_LENGTH;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.INV_SHORT_MAX_HOMOLOGY_HIGHER;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.INV_SHORT_MAX_HOMOLOGY_LOWER;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.INV_SHORT_MIN_AF_HIGHER;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.INV_SHORT_MIN_AF_LOWER;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.INV_SHORT_RATE_HIGHER;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.INV_SHORT_RATE_LOWER;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.MIN_AVG_FRAG_FACTOR;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.MIN_AVG_FRAG_STD_DEV_FACTOR;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.MIN_TRIMMED_ANCHOR_LENGTH;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.PON_INS_SEQ_FWD_STRAND_1;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.PON_INS_SEQ_FWD_STRAND_2;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.PON_INS_SEQ_REV_STRAND_1;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.PON_INS_SEQ_REV_STRAND_2;
import static com.hartwig.hmftools.esvee.common.FilterType.DUPLICATE;
import static com.hartwig.hmftools.esvee.common.FilterType.INV_SHORT_ISOLATED;
import static com.hartwig.hmftools.esvee.common.FilterType.MIN_ANCHOR_LENGTH;
import static com.hartwig.hmftools.esvee.common.FilterType.MIN_LENGTH;
import static com.hartwig.hmftools.esvee.common.FilterType.MIN_QUALITY;
import static com.hartwig.hmftools.esvee.common.FilterType.MIN_AF;
import static com.hartwig.hmftools.esvee.common.FilterType.MIN_SUPPORT;
import static com.hartwig.hmftools.esvee.common.FilterType.SGL;
import static com.hartwig.hmftools.esvee.common.FilterType.INV_SHORT_FRAG_LOW_VAF;
import static com.hartwig.hmftools.esvee.common.FilterType.DEL_SHORT_LOW_VAF;
import static com.hartwig.hmftools.esvee.common.FilterType.SHORT_FRAG_LENGTH;
import static com.hartwig.hmftools.esvee.common.FilterType.INV_SHORT_LOW_VAF_HOM;
import static com.hartwig.hmftools.esvee.common.FilterType.STRAND_BIAS;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sv.SvVcfTags;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;
import com.hartwig.hmftools.esvee.common.FilterType;
import com.hartwig.hmftools.esvee.common.FragmentLengthBounds;
import com.hartwig.hmftools.esvee.prep.types.DiscordantStats;

import htsjdk.variant.variantcontext.Genotype;

public class VariantFilters
{
    private final FilterConstants mFilterConstants;
    private final FragmentLengthBounds mFragmentLengthBounds;
    private final DiscordantStats mDiscordantStats;

    private final double mShortInversionRate;
    private final double mShortFragmentInversionRate;

    public VariantFilters(
            final FilterConstants filterConstants, final FragmentLengthBounds fragmentLengthBounds, final DiscordantStats discordantStats)
    {
        mFilterConstants = filterConstants;
        mFragmentLengthBounds = fragmentLengthBounds;
        mDiscordantStats = discordantStats;

        mShortInversionRate = mDiscordantStats.shortInversionRate();
        mShortFragmentInversionRate = mDiscordantStats.shortInversionFragmentRate();
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

        if(isInversionShortLowVafHomology(var))
            var.addFilter(INV_SHORT_LOW_VAF_HOM);

        if(isInversionShortFragmentLowVaf(var))
            var.addFilter(INV_SHORT_FRAG_LOW_VAF);

        if(isDeletionShortLowVaf(var))
            var.addFilter(DEL_SHORT_LOW_VAF);

        if(failsStrandBias(var))
            var.addFilter(STRAND_BIAS);
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

        return fragmentLevelBelowStatisticalLimit(var.averageFragmentLength(), var.splitFragmentCount());
    }

    private boolean fragmentLevelBelowStatisticalLimit(int svAvgLength, int splitFragCount)
    {
        if(svAvgLength == 0) // for now treat is this as a pass
            return false;

        int medianLength = mFragmentLengthBounds.Median;
        double stdDeviation = mFragmentLengthBounds.StdDeviation;

        double lowerStdDevThreshold = MIN_AVG_FRAG_STD_DEV_FACTOR * stdDeviation;
        double lowerLengthLimit = medianLength - max(MIN_AVG_FRAG_FACTOR * stdDeviation / sqrt(splitFragCount), lowerStdDevThreshold);

        return svAvgLength < lowerLengthLimit;
    }

    private boolean isInversionShortLowVafHomology(final Variant var)
    {
        // filter INVs if adjusted AF is low and length < 3K
        if(var.type() != INV)
            return false;

        if(var.adjustedLength() > INV_SHORT_LENGTH)
            return false;

        int inexactHomology = var.breakendStart().InexactHomology.length();

        if(inexactHomology < INV_SHORT_MAX_HOMOLOGY_LOWER)
            return false;

        double vafThreshold;

        if(inexactHomology >= INV_SHORT_MAX_HOMOLOGY_HIGHER)
        {
            vafThreshold = min(INV_SHORT_MIN_AF_HIGHER, mShortInversionRate * INV_SHORT_RATE_HIGHER);
        }
        else
        {
            vafThreshold = min(INV_SHORT_MIN_AF_LOWER, mShortInversionRate * INV_SHORT_RATE_LOWER);
        }

        return !anySamplesAboveAfThreshold(var, vafThreshold);
    }

    private boolean isInversionShortFragmentLowVaf(final Variant var)
    {
        // Filter for short inv (vaf < 0.05; length < 300; vaf/shortInvRate < 50)
        if(var.type() != INV)
            return false;

        if(var.adjustedLength() > INV_SHORT_FRAGMENT_LENGTH)
            return false;

        double vafThreshold = min(INV_SHORT_FRAGMENT_MIN_AF, INV_SHORT_FRAGMENT_AF_RATIO * mShortFragmentInversionRate);

        return !anySamplesAboveAfThreshold(var, vafThreshold);
    }

    private boolean isDeletionShortLowVaf(final Variant var)
    {
        // filter DELs if AF is < 0.05> and length < 3K
        if(var.type() != DEL)
            return false;

        if(var.adjustedLength() > DEL_ARTEFACT_SHORT_LENGTH)
            return false;

        int inexactHomology = var.breakendStart().InexactHomology.length();

        if(inexactHomology < DEL_ARTEFACT_MAX_HOMOLOGY)
            return false;

        int inferredFragmentLength = var.averageFragmentLength() + var.adjustedLength();

        double lengthThreshold = mFragmentLengthBounds.UpperBound * DEL_ARTEFACT_LENGTH_FACTOR;

        if(inferredFragmentLength >= lengthThreshold)
            return false;

        return !anySamplesAboveAfThreshold(var, DEL_ARTEFACT_MIN_AF);
    }

    private boolean failsStrandBias(final Variant var)
    {
        Breakend breakend = var.breakendStart();

        boolean hasPassing = false;
        boolean hasFailing = false;

        for(Genotype genotype : breakend.Context.getGenotypes())
        {
            int fragmentCount = breakend.fragmentCount(genotype);

            if(fragmentCount <= 1)
                continue;

            double strandBias = getGenotypeAttributeAsDouble(genotype, SvVcfTags.STRAND_BIAS, 0.5);

            String insSequence = var.insertSequence();

            if(strandBias == 1.0 && (insSequence.contains(PON_INS_SEQ_FWD_STRAND_1) || insSequence.contains(PON_INS_SEQ_FWD_STRAND_2)))
                hasFailing = true;
            else if(strandBias == 0 && (insSequence.contains(PON_INS_SEQ_REV_STRAND_1) || insSequence.contains(PON_INS_SEQ_REV_STRAND_2)))
                hasFailing = true;
            else
                hasPassing = true;
        }

        if(hasFailing && !hasPassing)
            return true;

        return false;
    }

    public void applyAdjacentFilters(final Map<String,List<Breakend>> chromosomeBreakends)
    {
        // filter INVs which are short, unchained and isolated from other variants (exlcuding other similar INVs)
        for(List<Breakend> breakends : chromosomeBreakends.values())
        {
            for(int i = 0; i < breakends.size(); ++i)
            {
                Breakend breakend = breakends.get(i);

                if(breakend.sv().filters().contains(INV_SHORT_ISOLATED))
                    continue;

                if(!isCandidateIsolatedInv(breakend.sv()))
                    continue;

                // search forwards and backwards for an adjacent non-artefact breakend
                boolean hasValidAdjacent = false;

                for(int j = 0; j <= 1; ++j)
                {
                    boolean searchUp = (j == 0);
                    int k = searchUp ? i + 1 : i - 1;

                    while(k >= 0 && k < breakends.size())
                    {
                        Breakend otherBreakend = breakends.get(k);

                        if(abs(otherBreakend.Position - breakend.Position) > INV_ADJACENT_LENGTH * 2) // stop the search
                            break;

                        if(hasAdjacentNonArtefact(breakend, otherBreakend))
                        {
                            hasValidAdjacent = true;
                            break;
                        }

                        k += searchUp ? 1 : -1;
                    }

                    if(hasValidAdjacent)
                        break;
                }

                if(hasValidAdjacent)
                    continue;

                breakend.sv().addFilter(INV_SHORT_ISOLATED);
            }
        }
    }

    private static boolean isCandidateIsolatedInv(final Variant var)
    {
        // a short INV which either isn't in a chain or has low UFP evidence
        return var.type() == INV && var.svLength() <= INV_ADJACENT_LENGTH
            && (!var.inChainedAssembly() || var.maxUniqueFragmentPositions() < INV_ADJACENT_MIN_UPS);
    }

    private static boolean hasAdjacentNonArtefact(final Breakend breakend, final Breakend otherBreakend)
    {
        if(otherBreakend == breakend.otherBreakend())
            return false;

        if(breakend.Orient == otherBreakend.Orient) // must be facing away
            return false;

        if(isCandidateIsolatedInv(otherBreakend.sv()))
            return false;

        // skip if the next breakend is also an artefact (including PON)
        if(otherBreakend.sv().filters().stream().anyMatch(x -> INV_ADJACENT_EXCLUDED_FILTERS.contains(x)))
            return false;

        // within close range
        int maxLowerPosition, minUpperPosition;
        if(breakend.Position < otherBreakend.Position)
        {
            maxLowerPosition = breakend.Position + breakend.InexactHomology.End;
            minUpperPosition = otherBreakend.Position + otherBreakend.InexactHomology.Start;
        }
        else
        {
            maxLowerPosition = otherBreakend.Position + otherBreakend.InexactHomology.End;
            minUpperPosition = breakend.Position + breakend.InexactHomology.Start;
        }

        return maxLowerPosition >= minUpperPosition - INV_ADJACENT_LENGTH;
    }

    public static void logFilterTypeCounts(final List<Variant> variantList)
    {
        Map<FilterType,Integer> filterCounts = Maps.newHashMap();
        int passingVariants = 0;

        for(Variant var : variantList)
        {
            if(var.filters().isEmpty())
            {
                ++passingVariants;
                continue;
            }

            for(FilterType filterType : var.filters())
            {
                int count = filterCounts.getOrDefault(filterType, 0);
                filterCounts.put(filterType, count + 1);
            }
        }

        SV_LOGGER.info("variants passing({}) filtered({})", passingVariants, variantList.size() - passingVariants);

        for(FilterType filterType : FilterType.values())
        {
            Integer count = filterCounts.get(filterType);

            if(count != null)
            {
                SV_LOGGER.debug("variant filter {}: count({})", filterType, count);
            }
        }
    }
}
