package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatPloidy;
import static com.hartwig.hmftools.linx.chaining.ChainPloidyLimits.ploidyOverlap;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.ADJACENT;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.CHAIN_SPLIT;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.FOLDBACK;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.NEAREST;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.PLOIDY_MATCH;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.PLOIDY_MAX;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.PLOIDY_OVERLAP;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.ONLY;
import static com.hartwig.hmftools.linx.chaining.ProposedLinks.PM_MATCHED;
import static com.hartwig.hmftools.linx.chaining.ProposedLinks.PM_NONE;
import static com.hartwig.hmftools.linx.chaining.ProposedLinks.PM_OVERLAP;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.hasLinkClash;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ChainRuleSelector
{
    private ChainFinder mChainFinder;

    // references from chain-finder
    private Map<SvBreakend, List<SvLinkedPair>> mSvBreakendPossibleLinks;
    private final Map<SvVarData, SvChainState> mSvConnectionsMap;
    private final List<SvLinkedPair> mSkippedPairs;
    private final List<SvVarData> mComplexDupCandidates;
    private final List<SvVarData> mFoldbacks;
    private final List<SvLinkedPair> mAdjacentMatchingPairs;
    private final List<SvLinkedPair> mAdjacentPairs;
    private final List<SvChain> mChains;

    private static final Logger LOGGER = LogManager.getLogger(ChainRuleSelector.class);

    public ChainRuleSelector(
            final ChainFinder chainFinder,
            final Map<SvBreakend, List<SvLinkedPair>> svBreakendPossibleLinks,
            final Map<SvVarData, SvChainState> svConnectionsMap,
            final List<SvLinkedPair> skippedPairs,
            final List<SvVarData> foldbacks,
            final List<SvVarData> complexDupCandidates,
            final List<SvLinkedPair> adjacentMatchingPairs,
            final List<SvLinkedPair> adjacentPairs,
            final List<SvChain> chains)
    {
        mChainFinder = chainFinder;
        mSvBreakendPossibleLinks = svBreakendPossibleLinks;
        mSvConnectionsMap = svConnectionsMap;
        mSkippedPairs = skippedPairs;
        mFoldbacks = foldbacks;
        mComplexDupCandidates = complexDupCandidates;
        mAdjacentMatchingPairs = adjacentMatchingPairs;
        mAdjacentPairs = adjacentPairs;
        mChains = chains;
    }

    public List<ProposedLinks> findSingleOptionPairs()
    {
        // find all breakends with only one other link options, or only to both breakends of the same SV - ie an INV
        List<ProposedLinks> proposedLinks = Lists.newArrayList();

        for(Map.Entry<SvBreakend, List<SvLinkedPair>> entry : mSvBreakendPossibleLinks.entrySet())
        {
            if(entry.getValue().isEmpty())
            {
                LOGGER.warn("breakend({}) has no possibles left, should be purged", entry.getKey().toString());
                continue;
            }

            if(entry.getValue().size() >= 2) // disable connections to an INV for now
                continue;

            SvBreakend limitingBreakend = entry.getKey();

            final SvLinkedPair newPair;

            if(entry.getValue().size() == 2)
            {
                // consider a link only to an INV as a single option
                final SvLinkedPair pair1 = entry.getValue().get(0);
                final SvLinkedPair pair2 = entry.getValue().get(1);
                final SvVarData otherSv1 = pair1.getOtherSV(limitingBreakend.getSV());
                final SvVarData otherSv2 = pair2.getOtherSV(limitingBreakend.getSV());

                if(otherSv1 == otherSv2 && otherSv1.type() == INV)
                {
                    newPair = pair1.length() < pair2.length() ? pair1 : pair2;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                newPair = entry.getValue().get(0);
            }

            // special case for DM DUPs - because they can link with themselves at the end, don't restrict their connectivity earlier on
            if(mChainFinder.isDoubleMinuteDup(limitingBreakend.getSV()))
                continue;

            if(mSkippedPairs.contains(newPair))
                continue;

            // avoid the second breakend also being a limiting factor and so adding this link twice
            if(proposedLinks.stream().map(x -> x.Links.get(0)).anyMatch(y -> y == newPair))
                continue;

            double ploidyFirst = mChainFinder.getUnlinkedBreakendCount(newPair.firstBreakend());
            double ploidySecond = mChainFinder.getUnlinkedBreakendCount(newPair.secondBreakend());

            if(ploidyFirst == 0 || ploidySecond == 0)
                continue;

            ProposedLinks proposedLink = new ProposedLinks(newPair, ONLY);
            proposedLink.addBreakendPloidies(newPair.firstBreakend(), ploidyFirst, newPair.secondBreakend(), ploidySecond);

            // check for another proposed link with a clashing breakend, and if found take the lower ploidy and short link
            boolean addNew = true;
            int index = 0;
            while(index < proposedLinks.size())
            {
                final ProposedLinks otherLink = proposedLinks.get(index);
                final SvLinkedPair otherPair = otherLink.Links.get(0);

                if(otherPair.hasLinkClash(newPair) || otherPair.oppositeMatch(newPair))
                {
                    if (copyNumbersEqual(otherLink.ploidy(), proposedLink.ploidy()))
                    {
                        if (proposedLink.shortestLinkDistance() > otherLink.shortestLinkDistance())
                        {
                            addNew = false;
                            break;
                        }

                        proposedLinks.remove(otherLink);
                        continue;
                    }
                    else if (proposedLink.ploidy() > otherLink.ploidy())
                    {
                        addNew = false;
                        break;
                    }
                    else
                    {
                        proposedLinks.remove(otherLink);
                        continue;
                    }
                }

                ++index;
            }

            if(addNew)
                proposedLinks.add(proposedLink);
        }

        return proposedLinks;
    }

    public List<ProposedLinks> findFoldbackPairs(final List<ProposedLinks> proposedLinks)
    {
        // both ends of a foldback connect to one end of another SV with ploidy >= 2x

        if(mFoldbacks.isEmpty())
            return proposedLinks;

        List<ProposedLinks> newProposedLinks = Lists.newArrayList();
        List<SvVarData> processedChainedFoldbacks = Lists.newArrayList(); // to avoid double-processing

        for(SvVarData foldback : mFoldbacks)
        {
            boolean isChainedFoldback = foldback.isChainedFoldback();

            if(isChainedFoldback)
            {
                if(processedChainedFoldbacks.contains(foldback))
                    continue;

                processedChainedFoldbacks.add(foldback);
                processedChainedFoldbacks.add(foldback.getChainedFoldbackSv());
            }

            SvBreakend foldbackStart = null;
            SvBreakend foldbackEnd = null;

            if(!isChainedFoldback)
            {
                foldbackStart = foldback.getBreakend(true);
                foldbackEnd = foldback.getBreakend(false);
            }
            else
            {
                if(foldback.getFoldbackBreakend(true) != null)
                {
                    foldbackStart = foldback.getBreakend(true);
                    foldbackEnd = foldback.getFoldbackBreakend(true);
                }
                else
                {
                    foldbackStart = foldback.getBreakend(false);
                    foldbackEnd = foldback.getFoldbackBreakend(false);
                }
            }

            int origFoldbackPloidy = foldback.getImpliedPloidy();

            double foldbackPloidy = min(
                    mChainFinder.getUnlinkedBreakendCount(foldbackStart),
                    mChainFinder.getUnlinkedBreakendCount(foldbackEnd));

            if(foldbackPloidy == 0)
                continue;

            List<SvLinkedPair> pairsOnFbStart = mSvBreakendPossibleLinks.get(foldbackStart);
            List<SvLinkedPair> pairsOnFbEnd = mSvBreakendPossibleLinks.get(foldbackEnd);

            if(pairsOnFbStart == null || pairsOnFbEnd == null)
                continue;

            // get the set of possible pairs for the start and end breakends of the foldback and look for:
            // a) the same non-foldback breakend linked to different breakends of the same foldback
            // b) at least a 2:1 ploidy ratio between the non-foldback breakend and the foldback breakend

            for(SvLinkedPair pairStart : pairsOnFbStart)
            {
                SvVarData nonFbVar = pairStart.getOtherSV(foldback);
                SvBreakend otherBreakend = pairStart.getOtherBreakend(foldbackStart);

                if(nonFbVar.getImpliedPloidy() < origFoldbackPloidy)
                    continue;

                double nonFoldbackPloidy = mChainFinder.getUnlinkedBreakendCount(otherBreakend);

                // foldback ploidy must be half or less to match for a potential chain split
                if(!copyNumbersEqual(foldbackPloidy * 2, nonFoldbackPloidy) && nonFoldbackPloidy <= foldbackPloidy)
                    continue;

                // does this exist in the other foldback breakend's set of possible pairs
                for(SvLinkedPair pairEnd : pairsOnFbEnd)
                {
                    SvBreakend otherBreakend2 = pairEnd.getOtherBreakend(foldbackEnd);

                    // check that available breakends support this SV being connected twice
                    if(otherBreakend != otherBreakend2)
                        continue;

                    SvChain targetChain = null;
                    for(SvChain chain : mChains)
                    {
                        if(!copyNumbersEqual(foldbackPloidy * 2, chain.ploidy()) && chain.ploidy() <= foldbackPloidy)
                            continue;

                        for(int se = SE_START; se <= SE_END; ++se)
                        {
                            boolean isStart = isStart(se);
                            SvBreakend chainBe = chain.getOpenBreakend(isStart);

                            if (chainBe == null)
                                continue;

                            if(chainBe == otherBreakend || chainBe == otherBreakend)
                            {
                                targetChain = chain;
                                break;
                            }
                        }
                    }

                    ProposedLinks proposedLink = new ProposedLinks(
                            Lists.newArrayList(pairStart, pairEnd), CHAIN_SPLIT, targetChain);

                    double linkedPloidy = targetChain != null ? targetChain.ploidy() : nonFoldbackPloidy;

                    proposedLink.addFoldbackBreakends(
                            foldbackStart, foldbackEnd, foldbackPloidy,
                            otherBreakend, nonFoldbackPloidy, linkedPloidy, otherBreakend.ploidyUncertainty());

                    newProposedLinks.add(proposedLink);

                    LOGGER.trace("foldback({}) ploidy({}) matched with breakend({}) ploidy({})",
                            foldback.id(), formatPloidy(foldbackPloidy), otherBreakend, formatPloidy(nonFoldbackPloidy));

                    break;
                }
            }
        }

        // now check for a match between any previously identified proposed links and this set
        // note that to preserve the additional proposed-link info for complex links, the higher rule is considered 'new'
        return restrictProposedLinks(newProposedLinks, proposedLinks, ONLY);
    }

    public List<ProposedLinks> findComplexDupPairs(List<ProposedLinks> proposedLinks)
    {
        // both ends of a foldback or complex DUP connect to one end of another SV with ploidy >= 2x
        if(mComplexDupCandidates.isEmpty())
            return proposedLinks;

        List<ProposedLinks> newProposedLinks = Lists.newArrayList();

        // the complex DUP SVs need to connect to both ends of either a single SV or a set of chained SVs with twice the ploidy
        // for a complex DUP of the form D - A - D, where A has double ploidy of D, and both ends of D connect to both ends of A
        for(SvVarData compDup : mComplexDupCandidates)
        {
            double compDupPloidy = mChainFinder.getUnlinkedCount(compDup);

            SvBreakend compDupBeStart = compDup.getBreakend(true);
            SvBreakend compDupBeEnd = compDup.getBreakend(false);

            List<SvLinkedPair> pairsOnCdStart = mSvBreakendPossibleLinks.get(compDupBeStart);
            List<SvLinkedPair> pairsOnCdEnd = mSvBreakendPossibleLinks.get(compDupBeEnd);

            if (pairsOnCdStart == null || pairsOnCdEnd == null)
                continue;

            // search existing chains for open chain ends match the set of possibles for the complex DUP and with twice the ploidy
            for(SvChain chain : mChains)
            {
                if(!copyNumbersEqual(compDupPloidy * 2, chain.ploidy()) && chain.ploidy() <= compDupPloidy)
                    continue;

                SvBreakend chainBeStart = chain.getOpenBreakend(true);
                SvBreakend chainBeEnd = chain.getOpenBreakend(false);

                if(chainBeStart == null || chainBeEnd == null)
                    continue;

                double chainStartPloidy = mChainFinder.getUnlinkedBreakendCount(chainBeStart);
                double chainEndPloidy = mChainFinder.getUnlinkedBreakendCount(chainBeEnd);

                if(chainStartPloidy == 0 || chainEndPloidy == 0)
                    continue;

                SvLinkedPair[] matchingPair = {null, null};

                for(int se = SE_START; se <= SE_END; ++se)
                {
                    boolean isStart = isStart(se);
                    SvBreakend chainBe = chain.getOpenBreakend(isStart);

                    if(chainBe == null)
                        continue;

                    // look for this link amongst the possible pairs
                    for(int se2 = SE_START; se2 <= SE_END; ++se2)
                    {
                        boolean isStart2 = isStart(se2);
                        List<SvLinkedPair> comDupPairs = isStart2 ? pairsOnCdStart : pairsOnCdEnd;
                        for (SvLinkedPair pair : comDupPairs)
                        {
                            if(pair.hasBreakend(chainBe))
                            {
                                matchingPair[se2] = pair;
                                break;
                            }
                        }
                    }
                }

                if(matchingPair[SE_START] == null || matchingPair[SE_END] == null)
                    continue;

                ProposedLinks proposedLink = new ProposedLinks(Lists.newArrayList(matchingPair[SE_START], matchingPair[SE_END]),
                        CHAIN_SPLIT, chain);

                proposedLink.addComDupBreakends(
                        compDupBeStart, compDupBeEnd, compDupPloidy,
                        chainBeStart, chainStartPloidy, chainBeEnd, chainEndPloidy, chain.ploidy(), chain.ploidyUncertainty());

                newProposedLinks.add(proposedLink);

                LOGGER.trace("comDup({}) ploidy({}) matched with chain breakends({} & {}) ploidy({})",
                        compDup.id(), formatPloidy(compDupPloidy),
                        chainBeStart.toString(), chainBeEnd.toString(), formatPloidy(chain.ploidy()));

                /*
                mDiagnostics.logCsv("COMP_DUP", compDup,
                        String.format("ploidy(%.1f-%.1f-%.1f) beStart(%s ploidy=%.1f-%.1f) beEnd(%s ploidy=%.1f-%.1f)",
                                compDup.ploidyMin(), compDup.ploidy(), compDup.ploidyMax(),
                                chainBeStart.toString(), chainBeStart.getSV().ploidyMin(), chainBeStart.getSV().ploidyMax(),
                                chainBeEnd.toString(), chainBeEnd.getSV().ploidyMin(), chainBeEnd.getSV().ploidyMax()));
                */
            }

            // alternatively check if the same SV can be linked to the complex DUP's breakends and satisifies the ploidy constraints
            for(SvLinkedPair pairStart : pairsOnCdStart)
            {
                SvVarData otherVar = pairStart.getOtherSV(compDup);
                SvBreakend otherBreakend = pairStart.getOtherBreakend(compDupBeStart);

                double otherBreakendPloidy = mChainFinder.getUnlinkedBreakendCount(otherBreakend);

                if(!copyNumbersEqual(compDupPloidy * 2, otherBreakendPloidy) && otherBreakendPloidy <= compDupPloidy)
                    continue;

                // does this exist in the other foldback breakend's set of possible pairs
                for(SvLinkedPair pairEnd : pairsOnCdEnd)
                {
                    SvBreakend otherBreakend2 = pairEnd.getOtherBreakend(compDupBeEnd);

                    // check that available breakends support this SV being connected twice
                    if(otherBreakend.getSV() != otherBreakend2.getSV())
                        continue;

                    double otherBreakendPloidy2 = mChainFinder.getUnlinkedBreakendCount(otherBreakend2);

                    if(copyNumbersEqual(compDupPloidy * 2, otherBreakendPloidy2) || otherBreakendPloidy2 > compDupPloidy)
                    {
                        ProposedLinks proposedLink = new ProposedLinks(Lists.newArrayList(pairStart, pairEnd),
                                CHAIN_SPLIT, null);

                        proposedLink.addComDupBreakends(
                                compDupBeStart, compDupBeEnd, compDupPloidy,
                                otherBreakend, otherBreakendPloidy, otherBreakend2, otherBreakendPloidy2,
                                otherBreakendPloidy, otherBreakend.ploidyUncertainty());

                        // proposedLink.addBreakendPloidies(compDupBeStart, compDupPloidy, otherBreakend, otherBreakendPloidy/2);
                        // proposedLink.addBreakendPloidies(compDupBeEnd, compDupPloidy, otherBreakend2, otherBreakendPloidy2/2, false);

                        newProposedLinks.add(proposedLink);

                        LOGGER.trace("comDup({}) ploidy({}) matched with breakends({} & {}) ploidy({} & {})",
                                compDup.id(), formatPloidy(compDupPloidy),
                                otherBreakend, otherBreakend2, formatPloidy(otherBreakendPloidy), formatPloidy(otherBreakendPloidy2));
                    }

                    break;
                }
            }
        }

        return restrictProposedLinks(newProposedLinks, proposedLinks, ONLY);
    }

    public List<ProposedLinks> findFoldbackToFoldbackPairs(final List<ProposedLinks> proposedLinks)
    {
        // look for 2 foldbacks facing each other
        if(mFoldbacks.isEmpty())
            return proposedLinks;

        List<ProposedLinks> newProposedLinks = Lists.newArrayList();

        // first gather up the foldback breakends as pairs
        List<SvBreakend> foldbackBreakends = Lists.newArrayList();

        for(SvVarData foldback : mFoldbacks)
        {
            if(foldback.getFoldbackBreakend(true) != null)
                foldbackBreakends.add(foldback.getBreakend(true));

            if(foldback.getFoldbackBreakend(false) != null)
                foldbackBreakends.add(foldback.getBreakend(false));
        }

        for(SvBreakend breakend : foldbackBreakends)
        {
            double foldbackPloidy = mChainFinder.getUnlinkedBreakendCount(breakend);

            if(foldbackPloidy == 0)
                continue;

            final List<SvLinkedPair> svLinks = mSvBreakendPossibleLinks.get(breakend);

            if(svLinks == null)
                continue;

            for(final SvLinkedPair pair : svLinks)
            {
                SvBreakend otherBreakend = pair.getOtherBreakend(breakend);

                if(!foldbackBreakends.contains(otherBreakend))
                    continue;

                double otherBreakendPloidy = mChainFinder.getUnlinkedBreakendCount(otherBreakend);

                if(mSkippedPairs.contains(pair))
                    continue;

                if(proposedLinks.stream().map(x -> x.Links.get(0)).anyMatch(y -> y == pair))
                    continue;

                LOGGER.trace("pair({}) of foldbacks with ploidy({} & {})",
                        pair.toString(), formatPloidy(foldbackPloidy), formatPloidy(otherBreakendPloidy));

                ProposedLinks proposedLink = new ProposedLinks(pair, FOLDBACK);
                proposedLink.addBreakendPloidies(breakend, foldbackPloidy, otherBreakend, otherBreakendPloidy);
                newProposedLinks.add(proposedLink);
            }
        }

        // now check for a match between any previously identified proposed links and this set
        return restrictProposedLinks(proposedLinks, newProposedLinks, FOLDBACK);
    }

    public List<ProposedLinks> findPloidyMatchPairs(List<ProposedLinks> proposedLinks)
    {
        List<ProposedLinks> newProposedLinks = Lists.newArrayList();

        if(!proposedLinks.isEmpty())
        {
            proposedLinks.stream().filter(x -> x.ploidyMatchType() == PM_MATCHED).forEach(x -> x.addRule(PLOIDY_MATCH));
            proposedLinks.stream().filter(x -> x.ploidyMatchType() == PM_OVERLAP).forEach(x -> x.addRule(PLOIDY_OVERLAP));
            return proposedLinks;
        }

        double currentMaxPloidy = 0;
        List<SvLinkedPair> addedLinks = Lists.newArrayList();

        for(SvChainState svConn : mSvConnectionsMap.values())
        {
            SvVarData var = svConn.SV;

            // check whether this SV has any possible links with SVs of the same (remaining) rep count
            for(int be = SE_START; be <= SE_END; ++be)
            {
                if(var.isNullBreakend() && be == SE_END)
                    continue;

                boolean isStart = isStart(be);
                double breakendPloidy = svConn.unlinked(be);

                if(breakendPloidy == 0)
                    continue;

                if(!copyNumbersEqual(breakendPloidy, currentMaxPloidy) && breakendPloidy < currentMaxPloidy)
                    continue;

                final SvBreakend breakend = var.getBreakend(isStart);
                final List<SvLinkedPair> svLinks = mSvBreakendPossibleLinks.get(breakend);

                if(svLinks == null)
                    continue;

                for(final SvLinkedPair pair : svLinks)
                {
                    if(addedLinks.contains(pair))
                        continue;

                    if(mSkippedPairs.contains(pair))
                        continue;

                    SvBreakend otherBreakend = pair.getOtherBreakend(breakend);

                    double otherBreakendPloidy = mChainFinder.getUnlinkedBreakendCount(otherBreakend);

                    String ploidyMatch = PM_NONE;
                    if(copyNumbersEqual(otherBreakendPloidy, breakendPloidy))
                    {
                        ploidyMatch = PM_MATCHED;
                    }
                    else if(ploidyOverlap(var.ploidyUncertainty(), breakendPloidy, otherBreakendPloidy, otherBreakend.ploidyUncertainty()))
                    {
                        ploidyMatch = PM_OVERLAP;
                    }
                    else
                    {
                        continue;
                    }

                    // LOGGER.trace("pair({}) with {} ploidy({} & {})",
                    //        pair.toString(), ploidyMatch, formatPloidy(breakendPloidy), formatPloidy(otherBreakendPloidy));

                    ProposedLinks proposedLink = (ploidyMatch == PM_MATCHED) ?
                            new ProposedLinks(pair, PLOIDY_MATCH) : new ProposedLinks(pair, PLOIDY_OVERLAP);

                    proposedLink.addBreakendPloidies(breakend, breakendPloidy, otherBreakend, otherBreakendPloidy);
                    newProposedLinks.add(proposedLink);

                    if(!copyNumbersEqual(proposedLink.ploidy(), currentMaxPloidy))
                        currentMaxPloidy = proposedLink.ploidy();
                }
            }
        }

        return newProposedLinks;
    }

    public List<ProposedLinks> findAdjacentMatchingPairs(List<ProposedLinks> proposedLinks)
    {
        if(mAdjacentMatchingPairs.isEmpty())
            return proposedLinks;

        if(!proposedLinks.isEmpty())
        {
            for(ProposedLinks proposedLink : proposedLinks)
            {
                // check for any links which are in the adjacent set
                if(proposedLink.Links.stream().anyMatch(x -> mAdjacentMatchingPairs.contains(x)))
                {
                    proposedLink.addRule(ADJACENT);
                }
            }

            return proposedLinks;
        }

        List<ProposedLinks> newProposedLinks = Lists.newArrayList();

        while(!mAdjacentMatchingPairs.isEmpty())
        {
            SvLinkedPair nextPair = mAdjacentMatchingPairs.get(0);

            mAdjacentMatchingPairs.remove(0);

            if(mChainFinder.matchesExistingPair(nextPair))
                continue;

            double ploidyFirst = mChainFinder.getUnlinkedBreakendCount(nextPair.firstBreakend());
            double ploidySecond = mChainFinder.getUnlinkedBreakendCount(nextPair.secondBreakend());

            if(ploidyFirst == 0 || ploidySecond == 0)
                continue;

            if(!copyNumbersEqual(ploidyFirst, ploidySecond))
                continue;

            // take the average ploidy or calculate a weighted ploidy already?
            // if these links have already been partially used, then incorrect to calculate a weighted ploidy
            ProposedLinks proposedLink = new ProposedLinks(nextPair, ADJACENT);
            proposedLink.addRule(PLOIDY_MATCH);
            proposedLink.addBreakendPloidies(nextPair.firstBreakend(), ploidyFirst, nextPair.secondBreakend(), ploidySecond);
            newProposedLinks.add(proposedLink);
        }

        return restrictProposedLinks(proposedLinks, newProposedLinks, ADJACENT);
    }

    public List<ProposedLinks> findAdjacentPairs(List<ProposedLinks> proposedLinks)
    {
        if(mAdjacentPairs.isEmpty())
            return proposedLinks;

        if(!proposedLinks.isEmpty())
        {
            for(ProposedLinks proposedLink : proposedLinks)
            {
                // check for any links which are in the adjacent set
                if(proposedLink.Links.stream().anyMatch(x -> mAdjacentPairs.contains(x)))
                {
                    proposedLink.addRule(ADJACENT);
                }
            }

            return proposedLinks;
        }

        List<ProposedLinks> newProposedLinks = Lists.newArrayList();

        while(!mAdjacentPairs.isEmpty())
        {
            SvLinkedPair nextPair = mAdjacentPairs.get(0);

            mAdjacentPairs.remove(0);

            if(mChainFinder.matchesExistingPair(nextPair))
                continue;

            double ploidyFirst = mChainFinder.getUnlinkedBreakendCount(nextPair.firstBreakend());
            double ploidySecond = mChainFinder.getUnlinkedBreakendCount(nextPair.secondBreakend());

            if(ploidyFirst == 0 || ploidySecond == 0)
                continue;

            // take the average ploidy or calculate a weighted ploidy already?
            // if these links have already been partially used, then incorrect to calculate a weighted ploidy

            ProposedLinks proposedLink = new ProposedLinks(nextPair, ADJACENT);
            proposedLink.addRule(PLOIDY_MATCH);
            proposedLink.addBreakendPloidies(nextPair.firstBreakend(), ploidyFirst, nextPair.secondBreakend(), ploidySecond);
            newProposedLinks.add(proposedLink);
        }

        return restrictProposedLinks(proposedLinks, newProposedLinks, ADJACENT);
    }

    public List<ProposedLinks> findHighestPloidy(List<ProposedLinks> proposedLinks)
    {
        List<ProposedLinks> newProposedLinks = Lists.newArrayList();

        if(!proposedLinks.isEmpty())
        {
            // take the highest from amongst the proposed links
            double maxPloidy = proposedLinks.stream().mapToDouble(x -> x.ploidy()).max().getAsDouble();

            proposedLinks.stream().filter(x -> copyNumbersEqual(maxPloidy, x.ploidy())).forEach(x -> x.addRule(PLOIDY_MAX));
            return proposedLinks;
        }

        double currentMaxPloidy = 0;
        List<SvLinkedPair> addedLinks = Lists.newArrayList();

        for(SvChainState svConn : mSvConnectionsMap.values())
        {
            SvVarData var = svConn.SV;

            // check whether this SV has any possible links with SVs of the same (remaining) rep count
            for(int be = SE_START; be <= SE_END; ++be)
            {
                if(var.isNullBreakend() && be == SE_END)
                    continue;

                boolean isStart = isStart(be);
                double breakendPloidy = svConn.unlinked(be);

                if(breakendPloidy == 0)
                    continue;

                if(!copyNumbersEqual(breakendPloidy, currentMaxPloidy) && breakendPloidy < currentMaxPloidy)
                    continue;

                final SvBreakend breakend = var.getBreakend(isStart);
                final List<SvLinkedPair> svLinks = mSvBreakendPossibleLinks.get(breakend);

                if(svLinks == null)
                    continue;

                for(final SvLinkedPair pair : svLinks)
                {
                    if(addedLinks.contains(pair))
                        continue;

                    if(mSkippedPairs.contains(pair))
                        continue;

                    SvBreakend otherBreakend = pair.getOtherBreakend(breakend);

                    double otherBreakendPloidy = mChainFinder.getUnlinkedBreakendCount(otherBreakend);

                    if(otherBreakendPloidy == 0)
                        continue;

                    double minPairPloidy = min(otherBreakendPloidy, breakendPloidy);

                    if(!copyNumbersEqual(currentMaxPloidy, minPairPloidy) && minPairPloidy < currentMaxPloidy)
                        continue;

                    currentMaxPloidy = max(minPairPloidy, currentMaxPloidy);

                    LOGGER.trace("pair({}) with max ploidy({} & {})",
                            pair.toString(), formatPloidy(breakendPloidy), formatPloidy(otherBreakendPloidy));

                    ProposedLinks proposedLink = new ProposedLinks(pair, PLOIDY_MAX);
                    proposedLink.addBreakendPloidies(breakend, breakendPloidy, otherBreakend, otherBreakendPloidy);
                    newProposedLinks.add(proposedLink);
                }
            }
        }

        return newProposedLinks;
    }

    public List<ProposedLinks> findNearest(List<ProposedLinks> proposedLinks)
    {
        // sorts links into shortest first and purges any conflicting longer links with conflicting breakends
        if(proposedLinks.isEmpty())
        {
            for (SvChainState svConn : mSvConnectionsMap.values())
            {
                SvVarData var = svConn.SV;

                // check whether this SV has any possible links with SVs of the same (remaining) rep count
                for (int be = SE_START; be <= SE_END; ++be)
                {
                    if (var.isNullBreakend() && be == SE_END)
                        continue;

                    boolean isStart = isStart(be);
                    double breakendPloidy = svConn.unlinked(be);

                    if (breakendPloidy == 0)
                        continue;

                    final SvBreakend breakend = var.getBreakend(isStart);
                    final List<SvLinkedPair> svLinks = mSvBreakendPossibleLinks.get(breakend);

                    if (svLinks == null)
                        continue;

                    for (final SvLinkedPair pair : svLinks)
                    {
                        if (mSkippedPairs.contains(pair))
                            continue;

                        SvBreakend otherBreakend = pair.getOtherBreakend(breakend);

                        double otherBreakendPloidy = mChainFinder.getUnlinkedBreakendCount(otherBreakend);

                        if (otherBreakendPloidy == 0)
                            continue;

                        ProposedLinks proposedLink = new ProposedLinks(pair, NEAREST);
                        proposedLink.addBreakendPloidies(breakend, breakendPloidy, otherBreakend, otherBreakendPloidy);
                        proposedLinks.add(proposedLink);
                    }
                }
            }
        }

        List<ProposedLinks> shortestLinks = Lists.newArrayList();

        for(final ProposedLinks proposedLink : proposedLinks)
        {
            int index = 0;
            boolean addNew = true;
            while(index < shortestLinks.size())
            {
                final ProposedLinks otherLink = shortestLinks.get(index);

                // for proposed links with any breakend clash, just keep the shortest
                if(hasLinkClash(proposedLink.Links, otherLink.Links))
                {
                    if(proposedLink.shortestLinkDistance() < otherLink.shortestLinkDistance())
                    {
                        shortestLinks.remove(index);
                    }
                    else
                    {
                        addNew = false;
                        break;
                    }
                }

                ++index;
            }

            if(!addNew)
                continue;

            // insert by shortest distance first
            index = 0;
            while(index < shortestLinks.size())
            {
                if(proposedLink.shortestLinkDistance() < shortestLinks.get(index).shortestLinkDistance())
                    break;

                ++index;
            }

            LOGGER.trace("adding shortest proposed link: {} index({})", proposedLink.toString(), index);

            shortestLinks.add(index, proposedLink);
        }

        if(shortestLinks.size() > 1)
        {
            LOGGER.trace("found {} shortest non-clashing proposed links", proposedLinks.size());
        }

        return shortestLinks;
    }

    private static List<ProposedLinks> restrictProposedLinks(
            List<ProposedLinks> proposedLinks, List<ProposedLinks> newProposedLinks, ChainingRule newRule)
    {
        // now check for a match between any previously identified proposed links and this set
        if(proposedLinks.isEmpty())
            return newProposedLinks;

        if(newProposedLinks.isEmpty())
            return proposedLinks;

        for(ProposedLinks proposedLink : proposedLinks)
        {
            for (ProposedLinks newProposedLink : newProposedLinks)
            {
                if(anyLinkMatch(proposedLink.Links, newProposedLink.Links))
                {
                    proposedLink.addRule(newRule);
                    break;
                }
            }
        }

        return proposedLinks;
    }

    private static boolean anyLinkMatch(final List<SvLinkedPair> links1, final List<SvLinkedPair> links2)
    {
        for(final SvLinkedPair pair : links1)
        {
            if(links2.contains(pair))
                return true;
        }

        return false;
    }

    public static void cullByPriority(List<ProposedLinks> proposedLinks)
    {
        if(proposedLinks.size() <= 1)
            return;

        // find the highest priority and remove any entries less than this
        int maxPriority = proposedLinks.stream().mapToInt(x -> x.priority()).max().getAsInt();

        int index = 0;
        while(index < proposedLinks.size())
        {
            if(proposedLinks.get(index).priority() < maxPriority)
                proposedLinks.remove(index);
            else
                ++index;
        }
    }

}
