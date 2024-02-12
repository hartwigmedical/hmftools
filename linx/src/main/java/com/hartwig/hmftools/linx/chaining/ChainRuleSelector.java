package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;
import static com.hartwig.hmftools.linx.chaining.ChainFinder.MIN_CHAINING_JCN_LEVEL;
import static com.hartwig.hmftools.linx.chaining.ChainJcnLimits.jcnMatch;
import static com.hartwig.hmftools.linx.chaining.ChainJcnLimits.jcnMatchForSplits;
import static com.hartwig.hmftools.linx.chaining.ChainJcnLimits.jcnOverlap;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.ADJACENT;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.ADJACENT_MATCH;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.CA_JCN_SUPPORT;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.COMP_DUP_SPLIT;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.FOLDBACK;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.FOLDBACK_SPLIT;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.NEAREST;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.JCN_MATCH;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.JCN_MAX;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.JCN_OVERLAP;
import static com.hartwig.hmftools.linx.chaining.ChainingRule.ONLY;
import static com.hartwig.hmftools.linx.chaining.FoldbackBreakendPair.addByPloidy;
import static com.hartwig.hmftools.linx.chaining.FoldbackBreakendPair.containsBreakendPair;
import static com.hartwig.hmftools.linx.chaining.FoldbackBreakendPair.removeBreakendPair;
import static com.hartwig.hmftools.linx.chaining.FoldbackBreakendPair.updateBreakendPair;
import static com.hartwig.hmftools.linx.chaining.ProposedLinks.PM_MATCHED;
import static com.hartwig.hmftools.linx.chaining.ProposedLinks.PM_NONE;
import static com.hartwig.hmftools.linx.chaining.ProposedLinks.PM_OVERLAP;
import static com.hartwig.hmftools.linx.types.LinkedPair.hasLinkClash;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.DbPair;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

public class ChainRuleSelector
{
    private ChainLinkAllocator mLinkAllocator;

    private int mClusterId;
    private boolean mHasReplication;
    private List<ChainingRule> mRulesToApply;

    private List<FoldbackBreakendPair> mFoldbackBreakendPairs; // a pair of breakends forming a breakend, either chained or single
    private boolean mFoldbacksInitialised;

    // references from chain-finder
    private final Map<SvBreakend,List<LinkedPair>> mSvBreakendPossibleLinks;
    private final ChainJcnLimits mJcnLimits;
    private final SvChainConnections mSvConnectionsMap;
    private final Map<SvVarData,List<LinkedPair>> mComplexDupCandidates;
    private final List<SvVarData> mFoldbacks;
    private final List<LinkedPair> mAdjacentMatchingPairs;
    private final List<LinkedPair> mAdjacentPairs;
    private final List<SvChain> mChains;

    public ChainRuleSelector(
            final ChainLinkAllocator linkAllocator,
            final ChainJcnLimits jcnLimits,
            final Map<SvBreakend, List<LinkedPair>> svBreakendPossibleLinks,
            final List<SvVarData> foldbacks,
            final Map<SvVarData,List<LinkedPair>> complexDupCandidates,
            final List<LinkedPair> adjacentMatchingPairs,
            final List<LinkedPair> adjacentPairs,
            final List<SvChain> chains)
    {
        mLinkAllocator = linkAllocator;
        mJcnLimits = jcnLimits;
        mSvBreakendPossibleLinks = svBreakendPossibleLinks;
        mSvConnectionsMap = linkAllocator.getSvConnections();
        mFoldbacks = foldbacks;
        mComplexDupCandidates = complexDupCandidates;
        mAdjacentMatchingPairs = adjacentMatchingPairs;
        mAdjacentPairs = adjacentPairs;
        mChains = chains;
        mRulesToApply = Lists.newArrayList();
        mFoldbackBreakendPairs = Lists.newArrayList();
        mFoldbacksInitialised = false;
    }

    public void initialise(int clusterId, boolean clusterHasReplication)
    {
        mHasReplication = clusterHasReplication;
        mClusterId = clusterId;

        mRulesToApply.clear();

        if(mHasReplication)
        {
            mRulesToApply.add(FOLDBACK_SPLIT);
            mRulesToApply.add(COMP_DUP_SPLIT);
            mRulesToApply.add(ONLY);

            mRulesToApply.add(JCN_MATCH);
            mRulesToApply.add(ADJACENT);
            mRulesToApply.add(JCN_MAX);
        }
        else
        {
            mRulesToApply.add(ONLY);
            mRulesToApply.add(ADJACENT_MATCH);
        }

        mRulesToApply.add(NEAREST);

        mFoldbacksInitialised = false;
        mFoldbackBreakendPairs.clear();
    }

    public List<ProposedLinks> findProposedLinks()
    {
        // find the next set of possible links to make according to the priority scheme
        // which is expressed in the set of chaining rules (ie the enumerated type)

        // some special cases
        List<ProposedLinks> proposedLinks = Lists.newArrayList();

        for(int i = 0; i < mRulesToApply.size(); ++i)
        {
            final ChainingRule rule = mRulesToApply.get(i);

            switch (rule)
            {
                case ONLY:
                    proposedLinks = findSingleOptionPairs(proposedLinks);
                    break;

                case FOLDBACK_SPLIT:
                    proposedLinks = findFoldbackLinks(proposedLinks);
                    break;

                case COMP_DUP_SPLIT:
                    proposedLinks = findComplexDups(proposedLinks);
                    break;

                case JCN_MATCH:
                    proposedLinks = findJcnMatchPairs(proposedLinks);
                    break;

                case ADJACENT:
                    proposedLinks = findAdjacentPairs(proposedLinks);
                    break;

                case ADJACENT_MATCH:
                    proposedLinks = findAdjacentMatchingPairs(proposedLinks);
                    break;

                case JCN_MAX:
                    proposedLinks = findHighestJcn(proposedLinks);
                    break;

                case NEAREST:
                    proposedLinks = findNearest(proposedLinks);
                    break;
            }

            if(proposedLinks.isEmpty())
                continue;

            mLinkAllocator.removeSkippedPairs(proposedLinks);

            // the top-priority rule acts like an annotation - whether the link violates the cluster ploidy across
            // one or more breakend segments - so even if a higher rule has only 1 proposed link, if it violates this condition
            // then keep searching for one which doesn't
            boolean linksHaveJcnSupport = checkClusterJcnSupport(proposedLinks);

            // the last rule 'nearest' can throw up multiple possible links but since they're not conflicting and ordered
            // from shortest to longest, there's no need to cull any more
            cullByPriority(proposedLinks);

            // if a single-top priority rule is left, also meeting the ploidy-support restriction, then take it
            if(proposedLinks.size() == 1 && linksHaveJcnSupport)
                return proposedLinks;
        }

        return proposedLinks;
    }

    public static class BreakendComparator implements Comparator<SvBreakend>
    {
        public int compare(final SvBreakend first, final SvBreakend second)
        {
            if(first.getSV().id() == second.getSV().id())
            {
                if(first.usesStart() == second.usesStart())
                    return 0;
                else
                    return first.usesStart() ? -1 : 1;
            }
            else
            {
                return first.getSV().id() < second.getSV().id() ? -1 : 1;
            }
        }
    }


    private List<ProposedLinks> findSingleOptionPairs(List<ProposedLinks> proposedLinks)
    {
        if(mSvBreakendPossibleLinks.isEmpty())
            return proposedLinks;

        // find all breakends with only one other link option
        if(!proposedLinks.isEmpty())
        {
            for(final ProposedLinks proposedLink : proposedLinks)
            {
                boolean hasSingleOption = false;

                for(final LinkedPair pair : proposedLink.Links)
                {
                    for(int se = SE_START; se <= SE_END; ++se)
                    {
                        List<LinkedPair> possiblePairs = mSvBreakendPossibleLinks.get(pair.getBreakend(isStart(se)));

                        if(possiblePairs == null)
                        {
                            LNX_LOGGER.error("breakend({}) has no possible pairs, from proposedLink: {}",
                                    pair.getBreakend(isStart(se)), proposedLink);
                            continue;
                        }

                        if(possiblePairs.size() == 1)
                        {
                            hasSingleOption = true;
                            break;
                        }
                    }

                    if(hasSingleOption)
                    {
                        proposedLink.addRule(ONLY);
                        break;
                    }
                }
            }
        }

        List<SvBreakend> sortedBreakends = mSvBreakendPossibleLinks.keySet().stream().sorted(new BreakendComparator()).collect(Collectors.toList());

        for(SvBreakend limitingBreakend : sortedBreakends)
        {
            List<LinkedPair> breakendPairs = mSvBreakendPossibleLinks.get(limitingBreakend);

            if(breakendPairs.isEmpty())
            {
                LNX_LOGGER.warn("breakend({}) has no possibles left, should be purged", limitingBreakend.toString());
                continue;
            }

            if(breakendPairs.size() >= 2) // disable connections to an INV for now
                continue;

            final LinkedPair newPair;

            if(breakendPairs.size() == 2)
            {
                // consider a link only to an INV as a single option
                final LinkedPair pair1 = breakendPairs.get(0);
                final LinkedPair pair2 = breakendPairs.get(1);
                final SvVarData otherSv1 = pair1.getOtherSV(limitingBreakend.getSV());
                final SvVarData otherSv2 = pair2.getOtherSV(limitingBreakend.getSV());

                if(otherSv1 == otherSv2 && otherSv1.type() == INV)
                {
                    newPair = pair1.positionDistance() < pair2.positionDistance() ? pair1 : pair2;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                newPair = breakendPairs.get(0);
            }

            if(mLinkAllocator.hasSkippedPairs(newPair))
                continue;

            // skip the duplicate link stored against the other breakend
            if(proposedLinks.stream().map(x -> x.Links.get(0)).anyMatch(y -> y == newPair))
                continue;

            double jcnFirst = mLinkAllocator.getUnlinkedBreakendCount(newPair.firstBreakend(), true);
            double jcnSecond = mLinkAllocator.getUnlinkedBreakendCount(newPair.secondBreakend(), true);

            if(jcnFirst == 0 || jcnSecond == 0)
                continue;

            ProposedLinks proposedLink = new ProposedLinks(newPair, ONLY);
            proposedLink.addBreakendPloidies(newPair.firstBreakend(), jcnFirst, newPair.secondBreakend(), jcnSecond);

            // check for another proposed link with a clashing breakend, and if found take the lower ploidy and short link
            boolean addNew = true;
            int index = 0;
            while(index < proposedLinks.size())
            {
                final ProposedLinks otherLink = proposedLinks.get(index);
                final LinkedPair otherPair = otherLink.Links.get(0);

                if(!otherPair.hasLinkClash(newPair) && !otherPair.oppositeMatch(newPair))
                {
                    ++index;
                    continue;
                }

                if(mHasReplication)
                {
                    if(copyNumbersEqual(otherLink.jcn(), proposedLink.jcn()))
                    {
                        if(proposedLink.shortestLinkDistance() > otherLink.shortestLinkDistance())
                        {
                            addNew = false;
                            break;
                        }

                        proposedLinks.remove(otherLink);
                    }

                    // keep both for now and let downstream rules decide
                    ++index;
                }
                else
                {
                    if(proposedLink.shortestLinkDistance() > otherLink.shortestLinkDistance())
                    {
                        addNew = false;
                        break;
                    }

                    proposedLinks.remove(otherLink);

                }
            }

            if(addNew)
                proposedLinks.add(proposedLink);
        }

        return proposedLinks;
    }

    private void updateFoldbackBreakends()
    {
        // start with previously identified foldbacks - from a single INV, a single-breakend foldback or a chained foldback
        // and then re-evaluate them as they are added to chains
        // eg if a foldback is used to split another SV, then the next chain itself continues to function as a foldback
        List<FoldbackBreakendPair> existingChainedPairs = Lists.newArrayList();

        if(!mFoldbacksInitialised)
        {
            // add non-chained simple INV foldbacks first - chained ones will be added next
            mFoldbacksInitialised = true;
            mFoldbackBreakendPairs.clear();

            for(SvVarData foldback : mFoldbacks)
            {
                if(foldback.isSingleBreakendFoldback() || foldback.isChainedFoldback())
                    continue;

                SvBreakend foldbackStart = foldback.getBreakend(true);
                SvBreakend foldbackEnd = foldback.getBreakend(false);

                double jcn = min(mLinkAllocator.getUnlinkedBreakendCount(foldbackStart), mLinkAllocator.getUnlinkedBreakendCount(foldbackEnd));

                if(jcn == 0)
                    continue;

                FoldbackBreakendPair fbPair = new FoldbackBreakendPair(foldbackStart, foldbackEnd, jcn, null);

                addByPloidy(mFoldbackBreakendPairs, fbPair);
            }
        }
        else
        {
            // remove any foldbacks with exhausted breakends or no more possible links
            int index = 0;
            while(index < mFoldbackBreakendPairs.size())
            {
                final FoldbackBreakendPair fbPair = mFoldbackBreakendPairs.get(index);

                final BreakendJcn bpStart = mLinkAllocator.getBreakendJcnData(fbPair.BreakendStart);
                final BreakendJcn bpEnd = mLinkAllocator.getBreakendJcnData(fbPair.BreakendEnd);

                double minUnlinkedPloidy = 0;
                boolean chainHasSimpleFoldback = false;

                if(fbPair.isChained())
                {
                    // if this foldback is still the open ends of a chain then it will have unexhausted ploidy and will be re-checked below
                    minUnlinkedPloidy = min(bpStart.unlinkedJcn(), bpEnd.unlinkedJcn());
                }
                else
                {
                    // a single-SV foldback must have unchained ploidy on both breakends
                    minUnlinkedPloidy = min(bpStart.UnchainedJcn, bpEnd.UnchainedJcn);

                    // and cannot have already been used as a foldback
                    if(!bpStart.Chains.isEmpty() && bpEnd.Chains.isEmpty())
                    {
                        chainHasSimpleFoldback = bpStart.Chains.stream().anyMatch(x -> chainHasFoldback(x, fbPair.BreakendStart.getSV()));
                    }
                }

                if(minUnlinkedPloidy < MIN_CHAINING_JCN_LEVEL || chainHasSimpleFoldback)
                {
                    LNX_LOGGER.debug("foldback pair({}) removed from consideration: minUnlinkedPloidy({}) chainHasSimpleFoldback({})",
                            fbPair, formatJcn(minUnlinkedPloidy), chainHasSimpleFoldback);
                    mFoldbackBreakendPairs.remove(index);
                    continue;
                }

                if(fbPair.isChained())
                    existingChainedPairs.add(fbPair); // will be checked against current chains next up
                else
                    fbPair.Ploidy = minUnlinkedPloidy;

                ++index;
            }
        }

        // also include chain ends which effectively form a foldback
        for(final SvChain chain : mChains)
        {
            final SvBreakend chainStart = chain.getOpenBreakend(true);
            final SvBreakend chainEnd = chain.getOpenBreakend(false);

            if(chainStart == null || chainEnd == null)
                continue;

            // are breakends consecutive and facing the same direction?
            if(!chainStart.chromosome().equals(chainEnd.chromosome()) || chainStart.orientation() != chainEnd.orientation())
                continue;

            // exclude if no possible pairs are permitted
            if(mSvBreakendPossibleLinks.get(chainStart) == null || mSvBreakendPossibleLinks.get(chainEnd) == null)
                continue;

            int startIndex = chainStart.getClusterChrPosIndex();
            int endIndex = chainEnd.getClusterChrPosIndex();

            if(chainStart != chainEnd && abs(startIndex - endIndex) != 1)
            {
                if(abs(startIndex - endIndex) > 2)
                    continue;

                // the foldback is invalid if it has a deletion bridge with overhang on the front-facing breakend
                boolean startIsFront = (chainStart.position() < chainEnd.position()) == (chainStart.orientation() == 1);
                final SvBreakend frontBreakend = startIsFront ? chainStart : chainEnd;

                // check for a DB which either allows or invalidates this foldback
                final DbPair dbLink = frontBreakend.getSV().getDBLink(frontBreakend.usesStart());

                if(dbLink != null && dbLink.length() > 0)
                    continue;
            }

            FoldbackBreakendPair fbPair = new FoldbackBreakendPair(chainStart, chainEnd, chain.jcn(), chain);

            boolean alreadyExists = containsBreakendPair(mFoldbackBreakendPairs, fbPair);

            removeBreakendPair(existingChainedPairs, fbPair);

            if(!alreadyExists)
            {
                LNX_LOGGER.debug("chain({}) adding chained foldback breakends({})", chain.id(), fbPair);
                addByPloidy(mFoldbackBreakendPairs, fbPair);
            }
            else
            {
                // update in case chain or ploidy has changed
                updateBreakendPair(mFoldbackBreakendPairs, fbPair);
            }
        }

        // remove any chained foldbacks no longer supported by the current set of chains
        existingChainedPairs.stream().forEach(x -> removeBreakendPair(mFoldbackBreakendPairs, x));
    }

    private static boolean chainHasFoldback(final SvChain chain, final SvVarData var)
    {
        for(int i = 0; i < chain.getLinkedPairs().size() - 1; ++i)
        {
            final LinkedPair pair = chain.getLinkedPairs().get(i);
            final LinkedPair nextPair = chain.getLinkedPairs().get(i + 1);

            if(pair.second() == var && nextPair.first() == var && pair.firstBreakend() == nextPair.secondBreakend())
                return true;
        }

        return false;
    }

    private static final int FOLDBACK_A_PRIORITY = 5;
    // private static final int FOLDBACK_B_FB_PRIORITY = 4;
    private static final int FOLDBACK_B_PRIORITY = 3;
    private static final int FOLDBACK_C_PRIORITY = 2;
    private static final int FOLDBACK_D_PRIORITY = 1;
    private static final int FOLDBACK_NO_PRIORITY = 0;

    private List<ProposedLinks> findFoldbackLinks(final List<ProposedLinks> proposedLinks)
    {
        updateFoldbackBreakends();

        if(mFoldbackBreakendPairs.isEmpty())
            return proposedLinks;

        // find the highest ploidy foldback where:
        //	a) it splits another chain with 2x ploidy
        //	b) it matches the ploidy of another chain
        //	c) it splits another foldback with >2x ploidy
        //	d) it is itself split by another foldback or complex duplication

        List<ProposedLinks> newProposedLinks = Lists.newArrayList();
        int linkScore = FOLDBACK_NO_PRIORITY; // keep track of the highest priority proposed links amongst the foldbacks

        double lastFoldbackPloidy = 0;

        // foldbacks are cached from highest to lowest ploidy already
        for(final FoldbackBreakendPair fbPair : mFoldbackBreakendPairs)
        {
            // break if the next foldback pair has a lower ploidy than any of those which have proposed a link
            if(!newProposedLinks.isEmpty() && !copyNumbersEqual(fbPair.Ploidy, lastFoldbackPloidy))
                break;

            lastFoldbackPloidy = fbPair.Ploidy;
            SvBreakend foldbackStart = fbPair.BreakendStart;
            SvBreakend foldbackEnd = fbPair.BreakendEnd;

            SvVarData foldback = foldbackStart.getSV();

            // if a foldback breakend is chained, use the ploidy of its chain
            SvChain foldbackChain = null;
            double foldbackPloidy = fbPair.Ploidy;
            double foldbackUncertainty = 0;

            if(fbPair.isChained())
            {
                foldbackChain = fbPair.Chain;
                foldbackUncertainty = foldbackChain.jcnUncertainty();
            }
            else
            {
                foldback.jcnUncertainty();
            }

            // search through all available pairs for the top priority type of link for this foldback
            List<LinkedPair> pairsOnFbStart = mSvBreakendPossibleLinks.get(foldbackStart);
            List<LinkedPair> pairsOnFbEnd = mSvBreakendPossibleLinks.get(foldbackEnd);

            if(pairsOnFbStart == null || pairsOnFbEnd == null)
                continue;

            pairsOnFbStart = Lists.newArrayList(pairsOnFbStart);
            pairsOnFbEnd = Lists.newArrayList(pairsOnFbEnd );

            cullDualOptionPairs(foldbackStart, pairsOnFbStart);

            for(LinkedPair pairStart : pairsOnFbStart)
            {
                SvVarData nonFbVar = pairStart.getOtherSV(foldback);
                SvBreakend otherBreakend = pairStart.getOtherBreakend(foldbackStart);

                // find the other pairing - would expect this to exist
                LinkedPair pairEnd = pairsOnFbEnd.stream()
                        .filter(x -> x.getOtherBreakend(foldbackEnd) == otherBreakend)
                        .findFirst().orElse(null);

                if(pairEnd == null)
                    continue;

                // first establish the available ploidy of this breakend and whether it's chained
                final BreakendJcn nonFbPloidyData = mLinkAllocator.getBreakendJcnData(otherBreakend);
                double nonFbPloidy = nonFbPloidyData.unlinkedJcn();

                // for low ploidy SVs, ploidy comparisons are imprecise, so don't allow splits
                boolean allowSplits = !(copyNumbersEqual(foldbackPloidy, nonFbPloidy) && foldbackPloidy < 1 && nonFbPloidy < 1);

                // a) first check if it splits another chain or SV with 2x ploidy
                if(allowSplits && jcnMatchForSplits(foldbackPloidy, foldbackUncertainty, nonFbPloidy, nonFbVar.jcnUncertainty())
                && !nonFbPloidyData.multiConnections())
                {
                    // a 2:1 splitting event
                    if(linkScore < FOLDBACK_A_PRIORITY)
                    {
                        // highest priority foldback type link, so clear any existing ones
                        newProposedLinks.clear();
                        linkScore = FOLDBACK_A_PRIORITY;
                    }

                    ProposedLinks proposedLink = new ProposedLinks(
                            Lists.newArrayList(pairStart, pairEnd), FOLDBACK_SPLIT, nonFbPloidyData.MaxJcnChain, foldbackChain);

                    proposedLink.addFoldbackBreakends(
                            foldbackStart, foldbackEnd, foldbackPloidy,
                            otherBreakend, nonFbPloidy, nonFbVar.jcnUncertainty());

                    LNX_LOGGER.trace("type-A: foldback breakends({} & {}) ploidy({}) exact split of breakend({}) ploidy({})",
                            foldbackStart, foldbackEnd, formatJcn(foldbackPloidy), otherBreakend, formatJcn(nonFbPloidy));

                    newProposedLinks.add(proposedLink);
                    continue;
                }

                if(linkScore >= FOLDBACK_A_PRIORITY)
                    continue;

                FoldbackBreakendPair otherFbPair = mFoldbackBreakendPairs.stream()
                        .filter(x -> x != fbPair)
                        .filter(x -> x.BreakendStart == otherBreakend || x.BreakendEnd == otherBreakend)
                        .findFirst().orElse(null);

                // b) check for an exact match with a chain or another SV
                if(!nonFbPloidyData.multiConnections()
                && jcnMatch(foldbackPloidy, foldbackUncertainty, nonFbPloidy, nonFbVar.jcnUncertainty()))
                {
                    if(linkScore < FOLDBACK_B_PRIORITY)
                    {
                        newProposedLinks.clear();
                        linkScore = FOLDBACK_B_PRIORITY;
                    }

                    // where there is a choice to be made between which breakends are used as is usually the case of foldbacks,
                    // choose the breakend with the most unlinked ploidy, which implies its other end is already chained
                    final ChainState fbConn = mLinkAllocator.getSvConnections().get(foldback);

                    final SvBreakend fbBreakend = fbConn.unlinked(foldbackStart.usesStart()) > fbConn.unlinked(foldbackEnd.usesStart())
                            ? foldbackStart : foldbackEnd;

                    ProposedLinks proposedLink = new ProposedLinks(LinkedPair.from(fbBreakend, otherBreakend), FOLDBACK);
                    proposedLink.addBreakendPloidies(fbBreakend, foldbackPloidy, otherBreakend, nonFbPloidy);

                    LNX_LOGGER.trace("type-B: foldback({}) ploidy({}) matched with {}({}) ploidy({})",
                            fbBreakend, formatJcn(foldbackPloidy), otherFbPair != null ? "foldback breakend" : "breakend",
                            otherBreakend, formatJcn(nonFbPloidy));

                    newProposedLinks.add(proposedLink);
                    continue;
                }

                if(linkScore >= FOLDBACK_B_PRIORITY)
                    continue;

                // c) check whether this foldback splits another foldback with > 2x ploidy

                // if this foldback splits another foldback then the other foldback's ploidy must be 2x or higher
                if(otherFbPair != null
                && (!jcnMatchForSplits(foldbackPloidy, foldbackUncertainty, nonFbPloidy, nonFbVar.jcnUncertainty())
                && nonFbPloidy > foldbackPloidy * 2))
                {
                    if(linkScore < FOLDBACK_C_PRIORITY) // 1
                    {
                        newProposedLinks.clear();
                        linkScore = FOLDBACK_C_PRIORITY; // 1
                    }

                    final SvChain targetChain =
                            nonFbPloidyData.MaxJcnChain != null && nonFbPloidy == nonFbPloidyData.MaxJcnChain.jcn() ?
                            nonFbPloidyData.MaxJcnChain : null;

                    ProposedLinks proposedLink = new ProposedLinks(
                            Lists.newArrayList(pairStart, pairEnd), FOLDBACK_SPLIT, targetChain, foldbackChain);

                    proposedLink.addFoldbackBreakends(
                            foldbackStart, foldbackEnd, foldbackPloidy,
                            otherBreakend, nonFbPloidy, nonFbVar.jcnUncertainty());

                    newProposedLinks.add(proposedLink);

                    LNX_LOGGER.trace("type-C: foldback breakends({} & {}) ploidy({}) non-exact split of foldback breakend({}) ploidy({})",
                            foldbackStart, foldbackEnd, formatJcn(foldbackPloidy), otherBreakend, formatJcn(nonFbPloidy));

                    continue;
                }

                // d) check whether the foldback is itself split by another foldback or complex duplication
                if(otherFbPair != null && allowSplits
                && jcnMatchForSplits(nonFbPloidy, nonFbVar.jcnUncertainty(), foldbackPloidy, foldbackUncertainty))
                {
                    linkScore = FOLDBACK_D_PRIORITY; // 0

                    // a 2:1 splitting event in reverse
                    ProposedLinks proposedLink = new ProposedLinks(
                            Lists.newArrayList(
                                    LinkedPair.from(otherFbPair.BreakendStart, foldbackStart),
                                    LinkedPair.from(otherFbPair.BreakendEnd, foldbackStart)),
                            FOLDBACK_SPLIT, foldbackChain, otherFbPair.Chain);

                    proposedLink.addFoldbackBreakends(
                            otherFbPair.BreakendStart, otherFbPair.BreakendEnd, otherFbPair.Ploidy,
                            foldbackStart, foldbackPloidy, foldbackUncertainty);

                    newProposedLinks.add(proposedLink);

                    LNX_LOGGER.trace("type-D: foldback breakend({}) ploidy({}) split by other foldback pair({})",
                            foldbackStart, formatJcn(foldbackPloidy), otherFbPair);

                    continue;
                }
            }
        }

            return newProposedLinks;
    }

    private void cullDualOptionPairs(final SvBreakend sourceBreakend, List<LinkedPair> pairs)
    {
        // if this breakend can connect to both ends of another SV (ie an INV), then removed the end which
        // is most exhausted (that is most chained)
        int i = 0;
        while(i < pairs.size())
        {
            final SvBreakend breakend = pairs.get(i).getOtherBreakend(sourceBreakend);

            boolean breakendRemoved = false;
            for(int j = i + 1; j < pairs.size(); ++j)
            {
                final SvBreakend breakend2 = pairs.get(j).getOtherBreakend(sourceBreakend);

                if(breakend.getSV() != breakend2.getSV())
                    continue;

                final ChainState svConn = mLinkAllocator.getSvConnections().get(breakend.getSV());

                if(svConn == null)
                    continue;

                if(copyNumbersEqual(svConn.unlinked(breakend.usesStart()), svConn.unlinked(breakend2.usesStart()))
                && pairs.get(i).positionDistance() > pairs.get(j).positionDistance())
                {
                    pairs.remove(i);
                    breakendRemoved = true;
                }
                else if(svConn.unlinked(breakend.usesStart()) < svConn.unlinked(breakend2.usesStart()))
                {
                    pairs.remove(i);
                    breakendRemoved = true;
                }
                else
                {
                    pairs.remove(j);
                }

                break;
            }

            if(!breakendRemoved)
                ++i;
        }
    }

    private List<ProposedLinks> findComplexDups(List<ProposedLinks> proposedLinks)
    {
        // check the complex dup candidates for any which can be proposed
        // for the ones around 2 different SVs, they must both be in the same chain and exhausted but for their ends
        if(mComplexDupCandidates.isEmpty())
            return proposedLinks;

        // only search for complex dups when no higher priority links have been found
        if(!proposedLinks.isEmpty())
            return proposedLinks;

        List<ProposedLinks> newProposedLinks = Lists.newArrayList();

        // the complex DUP SVs need to connect to both ends of either a single SV or a set of chained SVs with twice the ploidy
        // for a complex DUP of the form D - A - D, where A has double ploidy of D, and both ends of D connect to both ends of A
        List<SvVarData> invalidCompDups = Lists.newArrayList();

        for(Map.Entry<SvVarData,List<LinkedPair>> entry : mComplexDupCandidates.entrySet())
        {
            SvVarData compDup = entry.getKey();

            ChainState svConn = mSvConnectionsMap.get(compDup);
            if(svConn == null || svConn.hasConnections())
            {
                invalidCompDups.add(compDup);
                continue;
            }

            double compDupPloidy = svConn.Jcn;

            if(compDupPloidy == 0)
                continue;

            SvBreakend compDupBeStart = compDup.getBreakend(true);
            SvBreakend compDupBeEnd = compDup.getBreakend(false);

            List<LinkedPair> compDupLinks = entry.getValue();

            if(compDupLinks.size() != 2)
            {
                invalidCompDups.add(compDup);
                continue;
            }

            final LinkedPair pair1 = compDupLinks.get(0);
            final LinkedPair pair2 = compDupLinks.get(1);

            final SvVarData otherSv1 = pair1.getOtherSV(compDup);
            final SvVarData otherSv2 = pair2.getOtherSV(compDup);

            if(otherSv1 == otherSv2)
            {
                // a single SV duplicated by the complex DUP - check it isn't allocated yet
                ChainState otherSvConn = mSvConnectionsMap.get(otherSv1);

                final SvBreakend otherBreakend1 = otherSv1.getBreakend(true);
                final SvBreakend otherBreakend2 = otherSv1.getBreakend(false);

                if(otherSvConn == null || otherSvConn.hasConnections())
                {
                    invalidCompDups.add(compDup);
                    continue;
                }

                ProposedLinks proposedLink = new ProposedLinks(Lists.newArrayList(pair1, pair2), COMP_DUP_SPLIT, null, null);

                proposedLink.addComDupBreakends(
                        compDupBeStart, compDupBeEnd, compDupPloidy,
                        otherBreakend1, otherBreakend2, otherSvConn.Jcn);

                newProposedLinks.add(proposedLink);

                LNX_LOGGER.trace("comDup({}) ploidy({}) matched with single SV breakends({} & {}) ploidy({})",
                        compDup.id(), formatJcn(compDupPloidy), otherBreakend1,
                        otherBreakend2, formatJcn(otherSvConn.Jcn));
            }
            else
            {
                final SvBreakend otherBreakend1 = pair1.hasBreakend(compDupBeStart) ?
                        pair1.getOtherBreakend(compDupBeStart) : pair1.getOtherBreakend(compDupBeEnd);

                final SvBreakend otherBreakend2 = pair2.hasBreakend(compDupBeStart) ?
                        pair2.getOtherBreakend(compDupBeStart) : pair2.getOtherBreakend(compDupBeEnd);

                // search existing chains for open chain ends match the set of possibles for the complex DUP and with twice the ploidy
                for(SvChain chain : mChains)
                {
                    SvBreakend chainBeStart = chain.getOpenBreakend(true);
                    SvBreakend chainBeEnd = chain.getOpenBreakend(false);

                    if(!jcnMatchForSplits(compDupPloidy, compDup.jcnUncertainty(), chain.jcn(), chain.jcnUncertainty()))
                        continue;

                    if((chainBeStart == otherBreakend1 && chainBeEnd == otherBreakend2)
                    || (chainBeStart == otherBreakend2 && chainBeEnd == otherBreakend1))
                    {
                        ProposedLinks proposedLink = new ProposedLinks(Lists.newArrayList(pair1, pair2), COMP_DUP_SPLIT, chain, null);

                        proposedLink.addComDupBreakends(
                                compDupBeStart, compDupBeEnd, compDupPloidy,
                                chainBeStart, chainBeEnd, chain.jcn());

                        newProposedLinks.add(proposedLink);

                        LNX_LOGGER.trace("comDup({}) ploidy({}) matched with chain breakends({} & {}) ploidy({})",
                                compDup.id(), formatJcn(compDupPloidy),
                                chainBeStart.toString(), chainBeEnd.toString(), formatJcn(chain.jcn()));

                        break;
                    }
                }
            }

            /*
            mDiagnostics.logCsv("COMP_DUP", compDup,
                    String.format("ploidy(%.1f-%.1f-%.1f) beStart(%s ploidy=%.1f-%.1f) beEnd(%s ploidy=%.1f-%.1f)",
                            compDup.ploidyMin(), compDup.ploidy(), compDup.ploidyMax(),
                            chainBeStart.toString(), chainBeStart.getSV().ploidyMin(), chainBeStart.getSV().ploidyMax(),
                            chainBeEnd.toString(), chainBeEnd.getSV().ploidyMin(), chainBeEnd.getSV().ploidyMax()));
            */
        }

        invalidCompDups.stream().forEach(x -> mComplexDupCandidates.remove(x));

        return newProposedLinks;
    }

    private List<ProposedLinks> findJcnMatchPairs(List<ProposedLinks> proposedLinks)
    {
        // find pairs of matching ploidy breakends, taking the shortest where multiple exist
        List<ProposedLinks> newProposedLinks = Lists.newArrayList();

        boolean hasPloidySupportLinks = anyLinksHaveJcnSupport(proposedLinks);

        if(!proposedLinks.isEmpty())
        {
            proposedLinks.stream().filter(x -> x.jcnMatchType() == PM_MATCHED).forEach(x -> x.addRule(JCN_MATCH));
            proposedLinks.stream().filter(x -> x.jcnMatchType() == PM_OVERLAP).forEach(x -> x.addRule(JCN_OVERLAP));

            if(hasPloidySupportLinks)
                return proposedLinks;
        }

        List<LinkedPair> addedLinks = Lists.newArrayList();

        for(ChainState svConn : mSvConnectionsMap.values())
        {
            SvVarData var = svConn.SV;

            // check whether this SV has any possible links with SVs of the same (remaining) rep count
            for(int be = SE_START; be <= SE_END; ++be)
            {
                if(var.isSglBreakend() && be == SE_END)
                    continue;

                boolean isStart = isStart(be);
                double breakendPloidy = svConn.unlinked(be);

                if(breakendPloidy == 0)
                    continue;

                final SvBreakend breakend = var.getBreakend(isStart);
                final List<LinkedPair> svLinks = mSvBreakendPossibleLinks.get(breakend);

                if(svLinks == null)
                    continue;

                ProposedLinks bestProposedLink = null;

                for(final LinkedPair pair : svLinks)
                {
                    if(mLinkAllocator.hasSkippedPairs(pair))
                        continue;

                    if(addedLinks.contains(pair))
                        continue;

                    SvBreakend otherBreakend = pair.getOtherBreakend(breakend);

                    double otherBreakendPloidy = mLinkAllocator.getUnlinkedBreakendCount(otherBreakend, true);

                    if(otherBreakendPloidy == 0)
                        continue;

                    String ploidyMatch = PM_NONE;
                    if(copyNumbersEqual(otherBreakendPloidy, breakendPloidy))
                    {
                        ploidyMatch = PM_MATCHED;
                    }
                    else if(jcnOverlap(var.jcnUncertainty(), breakendPloidy, otherBreakendPloidy, otherBreakend.jcnUncertainty()))
                    {
                        ploidyMatch = PM_OVERLAP;
                    }
                    else
                    {
                        continue;
                    }

                    if(bestProposedLink != null && bestProposedLink.topRule() == JCN_OVERLAP && ploidyMatch != PM_MATCHED)
                    {
                        // no better and longer so keep searching for a better match
                        continue;
                    }

                    ProposedLinks proposedLink = (ploidyMatch == PM_MATCHED) ?
                            new ProposedLinks(pair, JCN_MATCH) : new ProposedLinks(pair, JCN_OVERLAP);

                    proposedLink.addBreakendPloidies(breakend, breakendPloidy, otherBreakend, otherBreakendPloidy);

                    bestProposedLink = proposedLink;
                    addedLinks.add(pair);

                    if(ploidyMatch == PM_MATCHED) // no need to keep searching through this breakend
                        break;
                }

                if(bestProposedLink != null)
                    newProposedLinks.add(bestProposedLink);
            }
        }

        checkClusterJcnSupport(newProposedLinks);

        // if a new link has ploidy support it will top anything found already (see earlier exit condition for proposed links)
        if(proposedLinks.isEmpty() || anyLinksHaveJcnSupport(newProposedLinks))
            return newProposedLinks;
        else
            return proposedLinks;
    }

    private List<ProposedLinks> findAdjacentMatchingPairs(List<ProposedLinks> proposedLinks)
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
            LinkedPair nextPair = mAdjacentMatchingPairs.get(0);

            mAdjacentMatchingPairs.remove(0);

            if(mLinkAllocator.matchesExistingPair(nextPair))
                continue;

            double ploidyFirst = mLinkAllocator.getUnlinkedBreakendCount(nextPair.firstBreakend(), true);
            double ploidySecond = mLinkAllocator.getUnlinkedBreakendCount(nextPair.secondBreakend(), true);

            if(ploidyFirst == 0 || ploidySecond == 0)
                continue;

            if(!copyNumbersEqual(ploidyFirst, ploidySecond))
                continue;

            // take the average ploidy or calculate a weighted ploidy already?
            // if these links have already been partially used, then incorrect to calculate a weighted ploidy
            ProposedLinks proposedLink = new ProposedLinks(nextPair, ADJACENT);
            proposedLink.addRule(JCN_MATCH);
            proposedLink.addBreakendPloidies(nextPair.firstBreakend(), ploidyFirst, nextPair.secondBreakend(), ploidySecond);
            newProposedLinks.add(proposedLink);
        }

        return restrictProposedLinks(proposedLinks, newProposedLinks, ADJACENT);
    }

    private List<ProposedLinks> findAdjacentPairs(List<ProposedLinks> proposedLinks)
    {
        if(mAdjacentPairs.isEmpty())
            return proposedLinks;

        boolean hasPloidySupportLinks = anyLinksHaveJcnSupport(proposedLinks);

        if(!proposedLinks.isEmpty())
        {
            // check for any links which are in the adjacent set
            proposedLinks.stream()
                    .filter(x -> x.Links.stream().anyMatch(y -> mAdjacentPairs.contains(y)))
                    .forEach(x -> x.addRule(ADJACENT));

            if(hasPloidySupportLinks)
                return proposedLinks;
        }

        List<ProposedLinks> newProposedLinks = Lists.newArrayList();

        while(!mAdjacentPairs.isEmpty())
        {
            LinkedPair nextPair = mAdjacentPairs.get(0);

            mAdjacentPairs.remove(0);

            if(mLinkAllocator.matchesExistingPair(nextPair))
                continue;

            double ploidyFirst = mLinkAllocator.getUnlinkedBreakendCount(nextPair.firstBreakend(), true);
            double ploidySecond = mLinkAllocator.getUnlinkedBreakendCount(nextPair.secondBreakend(), true);

            if(ploidyFirst == 0 || ploidySecond == 0)
                continue;

            ProposedLinks proposedLink = new ProposedLinks(nextPair, ADJACENT);
            proposedLink.addBreakendPloidies(nextPair.firstBreakend(), ploidyFirst, nextPair.secondBreakend(), ploidySecond);
            newProposedLinks.add(proposedLink);
        }

        checkClusterJcnSupport(newProposedLinks);

        // if a new link has ploidy support it will top anything found already (see earlier exit condition for proposed links)
        if(proposedLinks.isEmpty() || anyLinksHaveJcnSupport(newProposedLinks))
            return newProposedLinks;
        else
            return proposedLinks;

        // return restrictProposedLinks(proposedLinks, newProposedLinks, ADJACENT);
    }

    private List<ProposedLinks> findHighestJcn(List<ProposedLinks> proposedLinks)
    {
        List<ProposedLinks> newProposedLinks = Lists.newArrayList();

        boolean hasPloidySupportLinks = anyLinksHaveJcnSupport(proposedLinks);

        if(!proposedLinks.isEmpty())
        {
            // take the highest from amongst the proposed links
            double maxPloidy = proposedLinks.stream().mapToDouble(x -> x.jcn()).max().getAsDouble();

            proposedLinks.stream().filter(x -> copyNumbersEqual(maxPloidy, x.jcn())).forEach(x -> x.addRule(JCN_MAX));

            if(hasPloidySupportLinks)
                return proposedLinks;
        }

        double currentMaxPloidy = 0;
        List<LinkedPair> addedLinks = Lists.newArrayList();

        for(ChainState svConn : mSvConnectionsMap.values())
        {
            SvVarData var = svConn.SV;

            // check whether this SV has any possible links with SVs of the same (remaining) rep count
            for(int be = SE_START; be <= SE_END; ++be)
            {
                if(var.isSglBreakend() && be == SE_END)
                    continue;

                boolean isStart = isStart(be);
                double breakendPloidy = svConn.unlinked(be);

                if(breakendPloidy == 0)
                    continue;

                if(!copyNumbersEqual(breakendPloidy, currentMaxPloidy) && breakendPloidy < currentMaxPloidy)
                    continue;

                final SvBreakend breakend = var.getBreakend(isStart);
                final List<LinkedPair> svLinks = mSvBreakendPossibleLinks.get(breakend);

                breakendPloidy = mLinkAllocator.getUnlinkedBreakendCount(breakend, true);

                if(svLinks == null)
                    continue;

                for(final LinkedPair pair : svLinks)
                {
                    if(addedLinks.contains(pair))
                        continue;

                    if(mLinkAllocator.hasSkippedPairs(pair))
                        continue;

                    SvBreakend otherBreakend = pair.getOtherBreakend(breakend);

                    double otherBreakendPloidy = mLinkAllocator.getUnlinkedBreakendCount(otherBreakend, true);

                    if(otherBreakendPloidy == 0)
                        continue;

                    double minPairPloidy = min(otherBreakendPloidy, breakendPloidy);

                    if(!copyNumbersEqual(currentMaxPloidy, minPairPloidy) && minPairPloidy < currentMaxPloidy)
                        continue;

                    currentMaxPloidy = max(minPairPloidy, currentMaxPloidy);

                    LNX_LOGGER.trace("pair({}) with max ploidy({} & {})",
                            pair.toString(), formatJcn(breakendPloidy), formatJcn(otherBreakendPloidy));

                    ProposedLinks proposedLink = new ProposedLinks(pair, JCN_MAX);
                    proposedLink.addBreakendPloidies(breakend, breakendPloidy, otherBreakend, otherBreakendPloidy);
                    newProposedLinks.add(proposedLink);
                }
            }
        }

        checkClusterJcnSupport(newProposedLinks);

        // if a new link has ploidy support it will top anything found already (see earlier exit condition for proposed links)
        if(proposedLinks.isEmpty() || anyLinksHaveJcnSupport(newProposedLinks))
            return newProposedLinks;
        else
            return proposedLinks;
    }

    private List<ProposedLinks> findNearest(List<ProposedLinks> proposedLinks)
    {
        // sorts links into shortest first and purges any conflicting longer links with conflicting breakends
        if(proposedLinks.isEmpty())
        {
            for(ChainState svConn : mSvConnectionsMap.values())
            {
                SvVarData var = svConn.SV;

                // check whether this SV has any possible links with SVs of the same (remaining) rep count
                for(int be = SE_START; be <= SE_END; ++be)
                {
                    if(var.isSglBreakend() && be == SE_END)
                        continue;

                    boolean isStart = isStart(be);
                    double breakendPloidy = svConn.unlinked(be);

                    if(breakendPloidy == 0)
                        continue;

                    final SvBreakend breakend = var.getBreakend(isStart);

                    breakendPloidy = mLinkAllocator.getUnlinkedBreakendCount(breakend, true);

                    final List<LinkedPair> svLinks = mSvBreakendPossibleLinks.get(breakend);

                    if(svLinks == null)
                        continue;

                    for(final LinkedPair pair : svLinks)
                    {
                        if(mLinkAllocator.hasSkippedPairs(pair))
                            continue;

                        SvBreakend otherBreakend = pair.getOtherBreakend(breakend);

                        double otherBreakendPloidy = mLinkAllocator.getUnlinkedBreakendCount(otherBreakend, true);

                        if(otherBreakendPloidy == 0)
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
                else
                {
                    ++index;
                }
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

            LNX_LOGGER.trace("adding shortest proposed link: {} index({})", proposedLink.toString(), index);

            shortestLinks.add(index, proposedLink);
        }

        if(shortestLinks.size() > 1)
        {
            LNX_LOGGER.trace("found {} shortest non-clashing proposed links", shortestLinks.size());
        }

        return shortestLinks;
    }

    private boolean anyLinksHaveJcnSupport(final List<ProposedLinks> proposedLinks)
    {
        return proposedLinks.stream().anyMatch(x -> x.hasRule(CA_JCN_SUPPORT));
    }

    private boolean checkClusterJcnSupport(final List<ProposedLinks> proposedLinks)
    {
        boolean anyLinksHasClusterJcnSupport = false;
        for(ProposedLinks proposedLink : proposedLinks)
        {
            if(proposedLink.hasRule(CA_JCN_SUPPORT))
            {
                anyLinksHasClusterJcnSupport = true;
                continue;
            }

            if(proposedLink.Links.stream().anyMatch(x -> !mJcnLimits.linkHasJcnSupport(x, proposedLink.jcn())))
                continue;

            proposedLink.addRule(CA_JCN_SUPPORT);
            anyLinksHasClusterJcnSupport = true;
        }

        return anyLinksHasClusterJcnSupport;
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
            for(ProposedLinks newProposedLink : newProposedLinks)
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

    private static boolean anyLinkMatch(final List<LinkedPair> links1, final List<LinkedPair> links2)
    {
        for(final LinkedPair pair : links1)
        {
            if(links2.contains(pair))
                return true;
        }

        return false;
    }

    private void cullByPriority(List<ProposedLinks> proposedLinks)
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
