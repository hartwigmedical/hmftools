package com.hartwig.hmftools.esvee.vcfcompare.match;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.esvee.vcfcompare.CompareConfig;
import com.hartwig.hmftools.esvee.vcfcompare.common.VariantBreakend;

public class BreakendMatcher
{
    public final RefGenomeVersion mRefGenomeVersion;
    public final boolean mIncludeNonPass;

    private final List<BreakendMatch> mBreakendMatches = new ArrayList<>();

    public BreakendMatcher(RefGenomeVersion refGenomeVersion, boolean includeNonPass)
    {
        mRefGenomeVersion = refGenomeVersion;
        mIncludeNonPass = includeNonPass;
    }

    public BreakendMatcher(CompareConfig config)
    {
        mRefGenomeVersion = config.RefGenVersion;
        mIncludeNonPass = config.IncludeNonPass;
    }

    public void matchBreakends(
            Map<String, List<VariantBreakend>> oldChrBreakendMap,
            Map<String,List<VariantBreakend>> newChrBreakendMap,
            MatchType matchType,
            boolean checkOtherSide
    )
    {
        MatchFunctions.MatchFunction breakendMatcher = MatchType.getMatcher(matchType);

        int matchedCount = 0;

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mRefGenomeVersion.versionedChromosome(chromosome.toString());

            List<VariantBreakend> oldBreakends = oldChrBreakendMap.get(chrStr);
            List<VariantBreakend> newBreakends = newChrBreakendMap.get(chrStr);

            if(oldBreakends == null || newBreakends == null)
                continue;

            for(VariantBreakend oldBreakend : oldBreakends)
            {
                for(VariantBreakend newBreakend : newBreakends)
                {
                    if(oldBreakend.hasMatchedBreakend() || newBreakend.hasMatchedBreakend())
                        continue;

                    boolean hasMatch = breakendMatcher.match(oldBreakend, newBreakend, checkOtherSide);

                    if(hasMatch)
                    {
                        oldBreakend.MatchedBreakend = newBreakend;
                        newBreakend.MatchedBreakend = oldBreakend;

                        if(mIncludeNonPass || oldBreakend.isPassVariant() || newBreakend.isPassVariant())
                            mBreakendMatches.add(new BreakendMatch(oldBreakend, newBreakend, matchType));

                        matchedCount++;
                    }
                }
            }
        }

        if(matchedCount > 0)
        {
            SV_LOGGER.debug("Found {} variants with match type: {}", matchedCount, matchType);
        }
    }

    public void matchBreakends(
            Map<String,List<VariantBreakend>> oldChrBreakendMap,
            Map<String,List<VariantBreakend>> newChrBreakendMap,
            boolean checkOtherSide
    )
    {
        SV_LOGGER.info("Performing breakend matching using {} variants",
                mIncludeNonPass ? "ALL" : "PASS");

        matchBreakends(oldChrBreakendMap, newChrBreakendMap, MatchType.EXACT_MATCH, checkOtherSide);
        matchBreakends(oldChrBreakendMap, newChrBreakendMap, MatchType.COORDS_ONLY, checkOtherSide);
        matchBreakends(oldChrBreakendMap, newChrBreakendMap, MatchType.APPROX_MATCH, checkOtherSide);
    }

    public void matchBreakends(
            Map<String,List<VariantBreakend>> oldChrBreakendMap,
            Map<String,List<VariantBreakend>> newChrBreakendMap
    )
    {
        matchBreakends(oldChrBreakendMap, newChrBreakendMap, true);
    }

    private int gatherUnmatchedVariants(Map<String,List<VariantBreakend>> chrBreakendMap, boolean isOld)
    {
        int unmatchedVariantsCount = 0;

        for(List<VariantBreakend> breakends : chrBreakendMap.values())
        {
            if(breakends == null)
                continue;

            for(VariantBreakend breakend : breakends)
            {
                if(breakend.hasMatchedBreakend())
                    continue;

                if(mIncludeNonPass || breakend.isPassVariant())
                {
                    if(isOld)
                        mBreakendMatches.add(new BreakendMatch(breakend, null, MatchType.NO_MATCH));
                    else
                        mBreakendMatches.add(new BreakendMatch(null, breakend, MatchType.NO_MATCH));

                    unmatchedVariantsCount++;
                }
            }
        }

        return unmatchedVariantsCount;
    }

    public void gatherUnmatchedVariants(
            Map<String,List<VariantBreakend>> oldChrBreakendMap,
            Map<String,List<VariantBreakend>> newChrBreakendMap
    )
    {
        int unmatchedVariantsCount = 0;

        unmatchedVariantsCount += gatherUnmatchedVariants(oldChrBreakendMap, true);
        unmatchedVariantsCount += gatherUnmatchedVariants(newChrBreakendMap, false);

        if(unmatchedVariantsCount > 0)
        {
            SV_LOGGER.debug("Found {} unmatched variants", unmatchedVariantsCount);
        }
    }

    public List<BreakendMatch> getBreakendMatches(){ return mBreakendMatches; }
}
