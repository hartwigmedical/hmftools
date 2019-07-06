package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.linx.analysis.ChainPloidyLimits.CLUSTER_ALLELE_PLOIDY_MIN;
import static com.hartwig.hmftools.linx.analysis.ChainPloidyLimits.CLUSTER_AP;
import static com.hartwig.hmftools.linx.analysis.LinkFinder.getMinTemplatedInsertionLength;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LINK_TYPE_TI;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvChain;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/* ChainFinder - forms one or more chains from the SVs in a cluster

    Set-up:
    - form all assembled links and connect into chains
    - identify high-ploidy, foldbacks and complex DUP-type SVs, since these will be prioritised during chaining
    - create a cache of all possible linked pairs
    - create a cache of available breakends, including replicating each on accoriding to the SV's ploidy

    Replication count and ploidy
    - for clusters where all SVs have the same ploidy, none of the replication logic applies
    - in reality all SVs are integer multiples of each other according to their ploidy ratio, but in practice the ploidy min/max
    values need to be used to estimate replication

    Routine:
    - apply priority rules to find the next possible link(s)
    - add the link to an existing chain or a new chain if required
    - remove the breakends & link from further consideration
    - repeat until no further links can be made

    Optimisations:
    - for large clusters with many possible pairs on a chromosomal arm, only find the closest X initially (X = MaxPossiblePairs)
    - then when these possible pairs are exhausted for a breakend, search for more from the last pair's location
    - for foldbacks don't apply this restriction

    Priority rules:
    - Max-Replicated - find the SV(s) with the highest replication count, then select the one with the fewest possible links
    - Single-Option - if a breakend has only one possible link, select this one
    - Foldbacks - look for a breakend which can link to both ends of a foldback
    - Ploidy-Match - starting with the highest ploidy SV, only link SVs of the same ploidy
    - Resolving-SV - only link a high-ploidy SV to a lower one once all ploidy-match links are exhausted
    - Shortest - after all other rules, if there is more than 1 possible link then choose the shortest

    Rule selection
    1. Single-Option
    2. Foldbacks
    3. Ploidy-Match
    4. Max-Replicated (will possibly discard)
    5. Resolving-SV (possbly not relevant after rules 2 & 3 are exhausted)
    6. Shortest

*/

public class ChainFinder
{
    private int mClusterId;
    private boolean mHasReplication;

    // input state
    private List<SvVarData> mUniqueSVs;
    private List<SvVarData> mReplicatedSVs;
    private List<SvVarData> mFoldbacks;
    private boolean mIsClusterSubset;
    private List<SvLinkedPair> mAssembledLinks;
    private Map<String,List<SvBreakend>> mChrBreakendMap;

    // links a breakend to its position in a chromosome's breakend list - only applicable if a subset of SVs are being chained
    private Map<SvBreakend,Integer> mSubsetBreakendClusterIndexMap;

    // chaining state
    private List<SvLinkedPair> mSkippedPairs;
    private List<SvVarData> mComplexDupCandidates;
    private List<SvLinkedPair> mUniquePairs;
    private List<SvLinkedPair> mAdjacentMatchingPairs;
    private boolean mSkippedPair; // keep track of any excluded pair or SV without exiting the chaining routine

    private List<SvChain> mChains;
    private int mNextChainId;
    private List<SvVarData> mUnlinkedSVs;
    private Map<SvBreakend,List<SvBreakend>> mUnlinkedBreakendMap;
    private Map<SvVarData,Integer> mSvReplicationMap; // diminishes as SVs are added to chains
    private Map<SvVarData,Integer> mSvOriginalReplicationMap;
    private Map<SvBreakend,Integer> mBreakendLastLinkIndexMap;
    private ChainPloidyLimits mClusterPloidyLimits;

    private Map<SvBreakend,List<SvLinkedPair>> mSvBreakendPossibleLinks;

    private int mLinkIndex;
    private String mLinkReason;
    private boolean mIsValid;
    private boolean mLogVerbose;
    private boolean mRunValidation;
    private boolean mUseAllelePloidies;

    private static final Logger LOGGER = LogManager.getLogger(ChainFinder.class);

    // current unused
    private int mMaxPossibleLinks;
    private static int DEFAULT_MAX_POSSIBLE_LINKS = 0;

    public ChainFinder()
    {
        mUniqueSVs = Lists.newArrayList();
        mReplicatedSVs = Lists.newArrayList();
        mFoldbacks = Lists.newArrayList();
        mIsClusterSubset = false;
        mAssembledLinks = Lists.newArrayList();
        mChrBreakendMap = null;

        mClusterPloidyLimits = new ChainPloidyLimits();

        mAdjacentMatchingPairs = Lists.newArrayList();
        mSubsetBreakendClusterIndexMap = Maps.newHashMap();
        mBreakendLastLinkIndexMap = Maps.newHashMap();
        mComplexDupCandidates = Lists.newArrayList();
        mFoldbacks = Lists.newArrayList();
        mChains = Lists.newArrayList();
        mSkippedPairs = Lists.newArrayList();
        mSvBreakendPossibleLinks = Maps.newHashMap();
        mSvOriginalReplicationMap = Maps.newHashMap();
        mSvReplicationMap = Maps.newHashMap();
        mUnlinkedSVs = Lists.newArrayList();
        mUniquePairs = Lists.newArrayList();
        mUnlinkedBreakendMap = Maps.newHashMap();

        mHasReplication = false;
        mLogVerbose = false;
        mRunValidation = false;
        mIsValid = true;
        mSkippedPair = false;
        mNextChainId = 0;
        mLinkIndex = 0;
        mMaxPossibleLinks = DEFAULT_MAX_POSSIBLE_LINKS;
        mLinkReason = "";
        mUseAllelePloidies = false;
    }

    public void clear()
    {
        mClusterId = -1;
        mUniqueSVs.clear();
        mReplicatedSVs.clear();
        mFoldbacks.clear();
        mIsClusterSubset = false;
        mAssembledLinks.clear();
        mChrBreakendMap = null;

        mAdjacentMatchingPairs.clear();
        mSubsetBreakendClusterIndexMap.clear();
        mBreakendLastLinkIndexMap.clear();
        mComplexDupCandidates.clear();
        mFoldbacks.clear();
        mChains.clear();
        mSkippedPairs.clear();
        mSvBreakendPossibleLinks.clear();
        mSvOriginalReplicationMap.clear();
        mSvReplicationMap.clear();
        mUnlinkedBreakendMap.clear();
        mUniquePairs.clear();
        mUnlinkedSVs.clear();

        mNextChainId = 0;
        mLinkIndex = 0;
        mIsValid = true;
        mSkippedPair = false;
    }

    public void initialise(SvCluster cluster)
    {
        // critical that all state is cleared before the next run
        clear();

        mClusterId = cluster.id();
        mUniqueSVs.addAll(cluster.getSVs());
        mReplicatedSVs.addAll(cluster.getSVs(true));
        mFoldbacks.addAll(cluster.getFoldbacks());
        mAssembledLinks.addAll(cluster.getAssemblyLinkedPairs());
        mChrBreakendMap = cluster.getChrBreakendMap();
        mHasReplication = cluster.hasReplicatedSVs();
        mIsClusterSubset = false;
    }

    public void initialise(SvCluster cluster, final List<SvVarData> svList)
    {
        clear();

        mIsClusterSubset = true;
        mClusterId = cluster.id();

        mChrBreakendMap = Maps.newHashMap();

        for (final Map.Entry<String, List<SvBreakend>> entry : cluster.getChrBreakendMap().entrySet())
        {
            final List<SvBreakend> breakendList = Lists.newArrayList();

            for (final SvBreakend breakend : entry.getValue())
            {
                if (svList.contains(breakend.getSV()))
                {
                    mSubsetBreakendClusterIndexMap.put(breakend, breakendList.size());
                    breakendList.add(breakend);
                }
            }

            if (!breakendList.isEmpty())
            {
                mChrBreakendMap.put(entry.getKey(), breakendList);
            }
        }

        for(SvVarData var : svList)
        {
            mReplicatedSVs.add(var);

            if(var.isReplicatedSv())
            {
                mHasReplication = true;
                continue;
            }

            mUniqueSVs.add(var);

            if(var.isFoldback() && mFoldbacks.contains(var))
                mFoldbacks.add(var);

            for(int se = SE_START; se <= SE_END; ++se)
            {
                // only add an assembled link if it has a partner in the provided SV set, and can be replicated equally
                for (SvLinkedPair link : var.getAssembledLinkedPairs(isStart(se)))
                {
                    final SvVarData otherVar = link.getOtherSV(var);

                    if(!svList.contains(otherVar))
                        continue;

                    int maxRepCount = mHasReplication ? min(max(var.getReplicatedCount(),1), max(otherVar.getReplicatedCount(),1)) : 1;

                    long currentLinkCount = mAssembledLinks.stream().filter(x -> x.matches(link)).count();

                    if(currentLinkCount < maxRepCount)
                    {
                        mAssembledLinks.add(link);
                    }
                }
            }
        }
    }

    public void setLogVerbose(boolean toggle)
    {
        mLogVerbose = toggle;
        setRunValidation(toggle);
    }

    public void setRunValidation(boolean toggle) { mRunValidation = toggle; }
    public void setUseAllelePloidies(boolean toggle) { mUseAllelePloidies = toggle; }

    public final List<SvChain> getChains() { return mChains; }
    public double getValidAllelePloidySegmentPerc() { return mClusterPloidyLimits.getValidAllelePloidySegmentPerc(); }

    public void formClusterChains(boolean assembledLinksOnly)
    {
        if(mReplicatedSVs.size() < 2)
            return;

        if (mUniqueSVs.size() >= 4)
        {
            LOGGER.debug("cluster({}) starting chaining with assemblyLinks({}) svCount({} rep={})",
                    mClusterId, mAssembledLinks.size(), mUniqueSVs.size(), mReplicatedSVs.size());
        }

        mClusterPloidyLimits.initialise(mClusterId, mChrBreakendMap);

        buildChains(assembledLinksOnly);

        // report on how effective chaining was
        if(LOGGER.isDebugEnabled())
        {
            int breakendCount = (int)mUnlinkedBreakendMap.values().stream().count();

            List<SvVarData> uniqueUnlinkedSVs = Lists.newArrayList();

            for(final SvVarData var : mUnlinkedSVs)
            {
                if(!uniqueUnlinkedSVs.contains(var.getOrigSV()))
                    uniqueUnlinkedSVs.add(var.getOrigSV());
            }

            LOGGER.debug("cluster({}) chaining finished: chains({} links={}) unlinked SVs({} unique={}) breakends({} reps={}) validAllelePerc({})",
                    mClusterId, mChains.size(), mLinkIndex, mUnlinkedSVs.size(), uniqueUnlinkedSVs.size(),
                    mUnlinkedBreakendMap.size(), breakendCount,
                    String.format("%.2f", mClusterPloidyLimits.getValidAllelePloidySegmentPerc()));
        }

        if(!mIsValid)
        {
            LOGGER.warn("cluster({}) chain finding failed", mClusterId);
            return;
        }
    }

    public void addChains(SvCluster cluster)
    {
        // add these chains to the cluster, but skip any which are identical to existing ones,
        // which can happen for clusters with replicated SVs
        mChains.stream().forEach(chain -> checkAddNewChain(chain, cluster));

        for(int i = 0; i < cluster.getChains().size(); ++i)
        {
            final SvChain chain = cluster.getChains().get(i);

            if(LOGGER.isDebugEnabled())
            {
                LOGGER.debug("cluster({}) added chain({}) with {} linked pairs:", mClusterId, chain.id(), chain.getLinkCount());
                chain.logLinks();
            }

            chain.setId(i); // set after logging so can compare with logging during building
        }
    }

    private void checkAddNewChain(final SvChain newChain, SvCluster cluster)
    {
        if(!mHasReplication || cluster.getChains().isEmpty())
        {
            cluster.addChain(newChain, false);
            return;
        }

        // any identical chains will have their replicated SVs entirely removed
        for(final SvChain chain : cluster.getChains())
        {
            if(chain.identicalChain(newChain))
            {
                boolean allReplicatedSVs = newChain.getSvCount(false) == 0;

                LOGGER.debug("cluster({}) skipping duplicate chain({}) vs origChain({}) all replicated({})",
                        mClusterId, newChain.id(), chain.id(), allReplicatedSVs);

                chain.addToReplicationCount();

                // remove these replicated SVs as well as the replicated chain
                if(allReplicatedSVs)
                {
                    for (final SvVarData var : newChain.getSvList())
                    {
                        cluster.removeReplicatedSv(var);
                    }
                }

                return;
            }
        }

        cluster.addChain(newChain, false);
    }

    private void buildChains(boolean assembledLinksOnly)
    {
        setUnlinkedBreakends();

        setSvReplicationCounts();

        // first make chains out of any assembly links
        addAssemblyLinksToChains();

        if(assembledLinksOnly)
            return;

        if(mUseAllelePloidies && mHasReplication)
            mClusterPloidyLimits.determineBreakendPloidies();

        determinePossibleLinks();

        int iterationsWithoutNewLinks = 0; // protection against loops

        while (true)
        {
            mSkippedPair = false;
            int lastAddedIndex = mLinkIndex;

            List<SvLinkedPair> possiblePairs = findPossiblePairs();

            if(possiblePairs.isEmpty())
            {
                if(!mSkippedPair)
                    break;
            }
            else
            {
                processPossiblePairs(possiblePairs);
            }

            if(lastAddedIndex == mLinkIndex)
            {
                ++iterationsWithoutNewLinks;

                if (iterationsWithoutNewLinks > 5)
                {
                    LOGGER.error("cluster({}) 5 iterations without adding a new link", mClusterId);
                    mIsValid = false;
                    break;
                }
            }
            else
            {
                iterationsWithoutNewLinks = 0;
            }

            checkProgress();
        }
    }

    private static String PP_METHOD_ONLY = "ONLY";

    private List<SvLinkedPair> findPossiblePairs()
    {
        // proceed through the link-finding methods in priority
        List<SvLinkedPair> possiblePairs = findSingleOptionPairs();

        if(!possiblePairs.isEmpty())
        {
            mLinkReason = PP_METHOD_ONLY;
            return possiblePairs;
        }

        if (mHasReplication)
        {
            possiblePairs = findDuplicationPairs();

            if (!possiblePairs.isEmpty())
            {
                mLinkReason = "FB_DUP";
                return possiblePairs;
            }
        }

        possiblePairs = findAdjacentMatchingPairs();
        if (!possiblePairs.isEmpty())
        {
            mLinkReason = "ADJAC";
            return possiblePairs;
        }

        if(mHasReplication)
        {
            possiblePairs = findPloidyMatchPairs();

            if (!possiblePairs.isEmpty())
            {
                mLinkReason = "PL_MAT";
                return possiblePairs;
            }

            possiblePairs = findMaxReplicationPairs();

            if (!possiblePairs.isEmpty())
            {
                mLinkReason = "PL_MAX";
                return possiblePairs;
            }
        }

        // try all remaining
        List<SvBreakend> breakendList = mUnlinkedBreakendMap.keySet().stream().collect(Collectors.toList());
        possiblePairs = findFewestOptionPairs(breakendList, false);

        if(!possiblePairs.isEmpty())
        {
            mLinkReason = "SHORT";
            return possiblePairs;
        }

        return possiblePairs;
    }

    private List<SvLinkedPair> findPloidyMatchPairs()
    {
        List<SvLinkedPair> possiblePairs = Lists.newArrayList();

        // find the highest pair of SVs with matching ploidy
        int maxRepCount = 1;

        for(Map.Entry<SvVarData,Integer> entry : mSvReplicationMap.entrySet())
        {
            int repCount = entry.getValue();
            SvVarData var = entry.getKey();

            if(repCount <= 1)
                continue;

            if(repCount >= maxRepCount)
            {
                List<SvLinkedPair> newPairs = Lists.newArrayList();

                // check whether this SV has any possible links with SVs of the same (remaining) rep count
                for(int be = SE_START; be <= SE_END; ++be)
                {
                    if(var.isNullBreakend() && be == SE_END)
                        continue;

                    boolean isStart = isStart(be);
                    SvBreakend breakend = var.getBreakend(isStart);
                    List<SvLinkedPair> svLinks = mSvBreakendPossibleLinks.get(breakend);

                    if(svLinks == null)
                        continue;

                    for(SvLinkedPair pair : svLinks)
                    {
                        if(possiblePairs.contains(pair))
                            continue;

                        if(mSkippedPairs.contains(pair))
                            continue;

                        SvVarData otherVar = pair.getOtherSV(var);
                        Integer otherRepCount = mSvReplicationMap.get(otherVar);
                        if(otherRepCount != null && otherRepCount == repCount)
                        {
                            log(LOG_LEVEL_VERBOSE, String.format("pair(%s) with matching high-rep count(%s)", pair.toString(), repCount));
                            newPairs.add(pair);
                        }
                    }
                }

                if(newPairs.isEmpty())
                    continue;

                if(repCount > maxRepCount)
                {
                    maxRepCount = repCount;
                    possiblePairs.clear();
                }

                possiblePairs.addAll(newPairs);
            }
        }

        removeSkippedPairs(possiblePairs);

        return possiblePairs;
    }

    private List<SvLinkedPair> findDuplicationPairs()
    {
        // both ends of a foldback or complex DUP connect to one end of another SV with ploidy >= 2x
        List<SvLinkedPair> possiblePairs = Lists.newArrayList();

        // are there any foldbacks where both ends should or can only connect to a single other high-ploidy breakend
        List<SvVarData> replicatingSVs = Lists.newArrayList(mFoldbacks);
        replicatingSVs.addAll(mComplexDupCandidates);

        if(replicatingSVs.isEmpty())
            return possiblePairs;

        for(SvVarData var : replicatingSVs)
        {
            int varPloidy = var.getImpliedPloidy();

            // collect up possible pairs for this foldback
            List<SvLinkedPair> varPairs = Lists.newArrayList();

            for(int be = SE_START; be <= SE_END; ++be)
            {
                boolean isStart = isStart(be);

                if(var.getFoldbackBreakend(isStart) != null)
                {
                    SvBreakend breakend = var.getBreakend(isStart);
                    List<SvLinkedPair> bePairs = mSvBreakendPossibleLinks.get(breakend);
                    if(bePairs != null)
                    {
                        varPairs.addAll(bePairs);
                    }
                }
            }

            // then check if any of these run into a breakend of higher ploidy, and if so, take the nearest pair set
            // eg for A - B - A, where A has double ploidy of B, and both ends of B connect to same end of A
            // in the logic below, the SV in question is the B variant
            List<SvVarData> otherSVs = Lists.newArrayList();
            for(SvLinkedPair pair : varPairs)
            {
                SvVarData otherSV = pair.getOtherSV(var); // the SV with equal or higher ploidy
                int otherVarPloidy = otherSV.getImpliedPloidy();

                if(varPloidy > otherVarPloidy)
                    continue;

                if(otherSVs.contains(otherSV))
                    continue;

                otherSVs.add(otherSV);

                // at this point select the nearer of the breakend pair and link them both to one end of the higher-ploidy SV
                int endsMatchedOnOtherVarStart = 0;
                int endsMatchedOnOtherVarEnd = 0;
                long linkOnOtherVarStartLength = 0;
                long linkOnOtherVarEndLength = 0;
                int maxLinkOnOtherVarStart = 0;
                int maxLinkOnOtherVarEnd = 0;
                List<SvLinkedPair> startLinks = Lists.newArrayList();
                List<SvLinkedPair> endLinks = Lists.newArrayList();

                for(SvLinkedPair otherPair : varPairs)
                {
                    int maxPairCount = getMaxUnlinkedPairCount(otherPair);

                    if(maxPairCount == 0)
                        continue;

                    if(otherPair.hasBreakend(otherSV, true) && getUnlinkedBreakendCount(otherSV.getBreakend(true)) > 1
                            && (otherPair.hasBreakend(var, true) || otherPair.hasBreakend(var, false)))
                    {
                        // link can be made from start of this var to both ends of the other SV, and record the shortest pair length
                        ++endsMatchedOnOtherVarStart;
                        maxLinkOnOtherVarStart = max(maxLinkOnOtherVarStart, maxPairCount);
                        startLinks.add(otherPair);

                        if (linkOnOtherVarStartLength == 0 || otherPair.length() < linkOnOtherVarStartLength)
                            linkOnOtherVarStartLength = otherPair.length();
                    }
                    else if(otherPair.hasBreakend(otherSV, false) && getUnlinkedBreakendCount(otherSV.getBreakend(false)) > 1
                            && (otherPair.hasBreakend(var, true) || otherPair.hasBreakend(var, false)))
                    {
                        ++endsMatchedOnOtherVarEnd;
                        maxLinkOnOtherVarEnd = max(maxLinkOnOtherVarEnd, maxPairCount);
                        endLinks.add(otherPair);

                        if (linkOnOtherVarEndLength == 0 || otherPair.length() < linkOnOtherVarEndLength)
                            linkOnOtherVarEndLength = otherPair.length();
                    }
                }

                if(endsMatchedOnOtherVarStart < 2 && endsMatchedOnOtherVarEnd < 2)
                    continue;

                log(LOG_LEVEL_VERBOSE, String.format("SV(%s) ploidy(%d) foldback dual links: start(links=%d maxLinks=%d matched=%d) start(links=%d maxLinks=%d matched=%d)",
                        var.id(), varPloidy, startLinks.size(), maxLinkOnOtherVarStart, endsMatchedOnOtherVarStart,
                        endLinks.size(), maxLinkOnOtherVarEnd, endsMatchedOnOtherVarEnd));

                // check for conflicts between this pairing and any others involving the SVs, and if found
                // go with the shortest (could change this to the highest)
                if(!possiblePairs.isEmpty())
                {
                    boolean replacePossibles = false;
                    boolean hasClashes = false;

                    for(SvLinkedPair otherPair : possiblePairs)
                    {
                        if (otherPair.hasVariant(var) || otherPair.hasVariant(otherSV))
                        {
                            hasClashes = true;
                            int shorterNewPairs = 0;

                            if (endsMatchedOnOtherVarStart == 2)
                            {
                                shorterNewPairs = (int) startLinks.stream().filter(x -> x.length() < otherPair.length()).count();
                            }

                            if (endsMatchedOnOtherVarEnd == 2)
                            {
                                shorterNewPairs += (int) endLinks.stream().filter(x -> x.length() < otherPair.length()).count();
                            }

                            if (shorterNewPairs > 0)
                            {
                                replacePossibles = true;
                                break;
                            }
                        }
                    }

                    if(replacePossibles)
                    {
                        possiblePairs.clear();
                    }
                    else if(hasClashes)
                    {
                        continue;
                    }
                }

                if(endsMatchedOnOtherVarStart == 2 && endsMatchedOnOtherVarEnd == 2)
                {
                    // prioritise potential link count followed by length
                    if(maxLinkOnOtherVarStart > maxLinkOnOtherVarEnd || linkOnOtherVarStartLength < linkOnOtherVarEndLength)
                    {
                        possiblePairs.addAll(startLinks);
                    }
                    else
                    {
                        possiblePairs.addAll(endLinks);
                    }
                }
                else if(endsMatchedOnOtherVarStart == 2)
                {
                    possiblePairs.addAll(startLinks);
                }
                else if(endsMatchedOnOtherVarEnd == 2)
                {
                    possiblePairs.addAll(endLinks);
                }
            }
        }

        if(!possiblePairs.isEmpty())
        {
            removeSkippedPairs(possiblePairs);
        }

        return possiblePairs;
    }

    private List<SvLinkedPair> findSingleOptionPairs()
    {
        List<SvLinkedPair> restrictedPairs = Lists.newArrayList();

        for(Map.Entry<SvBreakend,List<SvLinkedPair>> entry : mSvBreakendPossibleLinks.entrySet())
        {
            if(entry.getValue().size() != 1)
                continue;

            SvBreakend limitingBreakend = entry.getKey();

            // confirm that more possible pairs can't be found
            if(mMaxPossibleLinks > 0)
            {
                if(addMorePossibleLinks(limitingBreakend, true))
                    continue;
            }

            SvLinkedPair newPair = entry.getValue().get(0);

            if(mSkippedPairs.contains(newPair))
                continue;

            int minLinkCount = getMaxUnlinkedPairCount(newPair);

            // add this if it doesn't clash, and if it does then take the highest ploidy first, following by shortest
            int index = 0;
            boolean canAdd = true;
            while(index < restrictedPairs.size())
            {
                SvLinkedPair otherPair = restrictedPairs.get(index);

                if(otherPair == newPair)
                {
                    canAdd = false;
                    break;
                }

                if(!otherPair.hasLinkClash(newPair) && !otherPair.oppositeMatch(newPair))
                {
                    ++index;
                    continue;
                }

                int otherMinLinkCount = getMaxUnlinkedPairCount(otherPair);

                if(minLinkCount < otherMinLinkCount || newPair.length() < otherPair.length()) // || ploidySumNew > ploidySumOther
                {
                    restrictedPairs.remove(index);
                }
                else
                {
                    canAdd = false;
                    break;
                }
            }

            if(canAdd)
            {
                log(LOG_LEVEL_VERBOSE, String.format("single-option pair(%s) limited by breakend(%s)",
                        newPair.toString(), limitingBreakend.toString()));
                restrictedPairs.add(newPair);
            }
        }

        return restrictedPairs;
    }

    private List<SvLinkedPair> findAdjacentMatchingPairs()
    {
        if(mAdjacentMatchingPairs.isEmpty())
            return mAdjacentMatchingPairs;

        // take next one
        List<SvLinkedPair> possiblePairs = Lists.newArrayList();

        while(!mAdjacentMatchingPairs.isEmpty())
        {
            SvLinkedPair nextPair = mAdjacentMatchingPairs.get(0);

            mAdjacentMatchingPairs.remove(0);

            if(matchesExistingPair(nextPair))
            {
                continue;
            }

            // check breakends are still available
            if(getUnlinkedBreakendCount(nextPair.getFirstBreakend()) == 0 || getUnlinkedBreakendCount(nextPair.getSecondBreakend()) == 0)
            {
                continue;
            }

            possiblePairs.add(nextPair);
            return possiblePairs;
        }

        return possiblePairs;
    }

    private List<SvLinkedPair> findMaxReplicationPairs()
    {
        // look at the remaining SVs with replication at least 2 and those with the highest
        // remaining replication count and then those with fewest options
        List<SvLinkedPair> possiblePairs = Lists.newArrayList();

        // first check if there are SVs with a higher replication count, and if so favour these first
        List<SvVarData> maxRepSVs = !mSvReplicationMap.isEmpty() ? getMaxReplicationSVs() : null;

        if(maxRepSVs == null || maxRepSVs.isEmpty())
            return possiblePairs;

        List<SvBreakend> breakendList = Lists.newArrayList();

        for (final SvBreakend breakend : mUnlinkedBreakendMap.keySet())
        {
            if (maxRepSVs.contains(breakend.getSV()))
            {
                breakendList.add(breakend);
            }
        }

        if(breakendList.isEmpty())
            return possiblePairs;

        if (mLogVerbose)
        {
            for (SvVarData var : maxRepSVs)
            {
                LOGGER.debug("restricted to rep SV: {} repCount({})", var.id(), mSvReplicationMap.get(var));
            }
        }

        // next take the pairings with the least alternatives
        possiblePairs = findFewestOptionPairs(breakendList, true);

        if (possiblePairs.isEmpty())
        {
            // these high-replication SVs yielded no possible links so remove them from consideration
            for (final SvVarData var : maxRepSVs)
            {
                log(LOG_LEVEL_VERBOSE, String.format("cluster(%s) removing high-replicated SV(%s %s)",
                        mClusterId, var.posId(), var.type()));

                mSvReplicationMap.remove(var);
            }

            // mSkippedPair = true;
        }

        return possiblePairs;
    }

    private void removeSkippedPairs(List<SvLinkedPair> possiblePairs)
    {
        // some pairs are temporarily unavailable for use (eg those which would close a chain)
        // to to avoid continually trying to add them, keep them out of consideration until a new links is added
        if(mSkippedPairs.isEmpty())
            return;

        int index = 0;
        while(index < possiblePairs.size())
        {
            SvLinkedPair pair = possiblePairs.get(index);

            if(mSkippedPairs.contains(pair))
                possiblePairs.remove(index);
            else
                ++index;
        }
    }

    private List<SvLinkedPair> findFewestOptionPairs(List<SvBreakend> breakendList, boolean isRestricted)
    {
        // of these pairs, do some have less alternatives links which could be made than others
        // eg if high-rep SVs are A and B, and possible links have been found A-C, A-D and B-E,
        // then count how many links could be made between C, D and E to other SVs, and then select the least restrictive
        // say A-C is the only link C can make and B-D is the only link that D can make, then would want to make them both
        int minPairCount = 0;
        List<SvLinkedPair> minLinkPairs = Lists.newArrayList();

        for(SvBreakend breakend : breakendList)
        {
            List<SvLinkedPair> possiblePairs = mSvBreakendPossibleLinks.get(breakend);

            if (possiblePairs == null || possiblePairs.isEmpty())
                continue;

            if (minPairCount == 0 || possiblePairs.size() < minPairCount)
            {
                minLinkPairs.clear();
                minPairCount = possiblePairs.size();
            }

            if (possiblePairs.size() == minPairCount)
            {
                for (SvLinkedPair pair : possiblePairs)
                {
                    if (!minLinkPairs.contains(pair))
                        minLinkPairs.add(pair);
                }
            }

            if (isRestricted)
            {
                // also check this max-rep SV's pairings to test how they are restricted
                for (SvLinkedPair pair : possiblePairs)
                {
                    SvBreakend otherBreakend = pair.getOtherBreakend(breakend);

                    List<SvLinkedPair> otherPossiblePairs = mSvBreakendPossibleLinks.get(otherBreakend);

                    if (otherPossiblePairs == null || otherPossiblePairs.isEmpty())
                        continue; // logical assert

                    if (minPairCount == 0 || otherPossiblePairs.size() < minPairCount)
                    {
                        minLinkPairs.clear();
                        minPairCount = otherPossiblePairs.size();
                    }

                    if (otherPossiblePairs.size() == minPairCount)
                    {
                        for (SvLinkedPair otherPair : otherPossiblePairs)
                        {
                            if (!minLinkPairs.contains(otherPair))
                                minLinkPairs.add(otherPair);
                        }
                    }
                }
            }
        }

        removeSkippedPairs(minLinkPairs);

        return minLinkPairs;
    }

    private void processPossiblePairs(List<SvLinkedPair> possiblePairs)
    {
        // now the top candidates to link have been found, take the shortest of them and add this to a chain
        // where possible, add links multiple times according to the min replication of the breakends involved
        // after each link is added, check whether any breakend now has only one link option
        boolean linkAdded = false;

        boolean isRestrictedSet = mLinkReason == PP_METHOD_ONLY;

        while (!possiblePairs.isEmpty())
        {
            SvLinkedPair shortestPair = null;
            for (SvLinkedPair pair : possiblePairs)
            {
                log(LOG_LEVEL_VERBOSE, String.format("method(%s) possible pair: %s length(%s)",
                        mLinkReason, pair.toString(), pair.length()));

                if (shortestPair == null || pair.length() < shortestPair.length())
                {
                    shortestPair = pair;
                }
            }

            possiblePairs.remove(shortestPair);

            // log(LOG_LEVEL_VERBOSE, String.format("shortest possible pair: %s length(%s)", shortestPair.toString(), shortestPair.length()));

            int pairRepeatCount = 1;

            if(mHasReplication)
            {
                int beStartCount = getUnlinkedBreakendCount(shortestPair.getBreakend(true));
                int beEndCount = getUnlinkedBreakendCount(shortestPair.getBreakend(false));

                if(beStartCount > 1 && beEndCount > 1)
                {
                    pairRepeatCount = min(beStartCount, beEndCount);
                }

                if(pairRepeatCount > 1)
                {
                    LOGGER.debug("repeating pair({}) {} times", shortestPair.toString(), pairRepeatCount);
                }
            }

            for(int i = 0; i < pairRepeatCount; ++i)
            {
                linkAdded |= addPairToChain(shortestPair);
            }

            if(!mIsValid)
                return;

            // check whether after adding a link, some SV breakends have only a single possible link
            if(!isRestrictedSet)
            {
                List<SvLinkedPair> restrictedPairs = findSingleOptionPairs();
                if(!restrictedPairs.isEmpty())
                {
                    possiblePairs = restrictedPairs;
                    isRestrictedSet = true;
                }
            }

            if(!isRestrictedSet)
            {
                // having added a new pair, remove any other conflicting pairs
                int index = 0;
                while (index < possiblePairs.size())
                {
                    SvLinkedPair pair = possiblePairs.get(index);
                    if (pair.oppositeMatch(shortestPair))
                    {
                        possiblePairs.remove(index);
                        continue;
                    }

                    if (mHasReplication)
                    {
                        if (findUnlinkedMatchingBreakend(pair.getBreakend(true)) == null
                        || findUnlinkedMatchingBreakend(pair.getBreakend(false)) == null)
                        {
                            // replicated instances exhausted
                            possiblePairs.remove(index);
                            continue;
                        }
                    }
                    else if (pair.hasLinkClash(shortestPair))
                    {
                        possiblePairs.remove(index);
                        continue;
                    }

                    ++index;
                }
            }
        }

        if(linkAdded)
        {
            mSkippedPairs.clear(); // any skipped links can now be re-evaluated
        }
    }

    private static int SPEC_LINK_INDEX = -1;
    // private static int SPEC_LINK_INDEX = 26;

    private boolean addPairToChain(final SvLinkedPair pair)
    {
        if(mLinkIndex == SPEC_LINK_INDEX)
        {
            LOGGER.debug("specific link index({}) pair({})", mLinkIndex, pair.toString());
        }

        // attempt to add to existing chain
        boolean addedToChain = false;
        boolean[] pairToChain = {false, false};

        SvBreakend unlinkedBeFirst = null;
        SvBreakend unlinkedBeSecond = null;
        final SvLinkedPair newPair;

        if(!mHasReplication)
        {
            newPair = pair;
        }
        else
        {
            // this pair was created from the set of possibles, but needs to make use of unlinked breakends
            unlinkedBeFirst = findUnlinkedMatchingBreakend(pair.getBreakend(true));
            unlinkedBeSecond = findUnlinkedMatchingBreakend(pair.getBreakend(false));

            if(unlinkedBeFirst == null || unlinkedBeSecond == null)
            {
                // tolerate missed assembly links while a more robust approach is determined for ploidy discrepancies
                if(pair.isAssembled())
                {
                    LOGGER.warn("cluster({}) missed assembly link", mClusterId);
                    return false;
                }

                mIsValid = false;
                LOGGER.error("new pair breakendStart({} valid={}) and breakendEnd({} valid={}) no unlinked match found",
                        pair.getBreakend(true).toString(), unlinkedBeFirst != null,
                        pair.getBreakend(false).toString(), unlinkedBeSecond != null);

                return false;
            }

            newPair = new SvLinkedPair(unlinkedBeFirst.getSV(), unlinkedBeSecond.getSV(), LINK_TYPE_TI,
                    unlinkedBeFirst.usesStart(), unlinkedBeSecond.usesStart());

            if(pair.isAssembled())
                newPair.setIsAssembled();
        }

        boolean linkClosesChain = false;

        for(SvChain chain : mChains)
        {
            // test this link against each ends to the chain
            boolean addToStart = false;
            pairToChain[0] = pairToChain[1] = false; // reset for scenario where skipped adding to both ends of chain

            for(int be = SE_START; be <= SE_END; ++be)
            {
                boolean isStart = isStart(be);
                final SvVarData chainSV = chain.getChainEndSV(isStart);

                if (chain.canAddLinkedPair(newPair, isStart, true))
                {
                    addToStart = isStart;

                    if (chainSV.equals(newPair.first(), true))
                    {
                        pairToChain[SE_START] = true;

                        if(chainSV != newPair.first())
                        {
                            // the correct SV was matched, but a different instance, so switch it for one matching the chain end
                            newPair.replaceFirst(chainSV);
                        }
                    }
                    else
                    {
                        pairToChain[SE_END] = true;

                        if (chainSV != newPair.second())
                        {
                            newPair.replaceSecond(chainSV);
                        }
                    }
                }
            }

            if(!pairToChain[SE_START] && !pairToChain[SE_END])
            {
                continue;
            }

            if(pairToChain[SE_START] && pairToChain[SE_END])
            {
                // the link can be added to both ends, which would close the chain - so search for an alternative SV on either end
                // to keep it open while still adding the link
                boolean replacementFound = false;

                if(mHasReplication)
                {
                    for (int be = SE_START; be <= SE_END; ++be)
                    {
                        boolean isStart = isStart(be);

                        SvBreakend openBreakend = chain.getOpenBreakend(isStart); // this will be one of the pair breakends

                        if(openBreakend == null)
                            continue; // eg ending on a SGL

                        List<SvBreakend> possibleBreakends = mUnlinkedBreakendMap.get(openBreakend.getOrigBreakend());
                        SvVarData chainSV = chain.getChainEndSV(isStart);

                        if (possibleBreakends == null || possibleBreakends.isEmpty())
                            continue;

                        for (SvBreakend otherBreakend : possibleBreakends)
                        {
                            if (otherBreakend.getSV() != chainSV)
                            {
                                replacementFound = true;

                                if (newPair.first() == chainSV)
                                    newPair.replaceFirst(otherBreakend.getSV());
                                else
                                    newPair.replaceSecond(otherBreakend.getSV());

                                pairToChain[be] = false;
                                addToStart = !isStart;
                                break;
                            }
                        }

                        if (replacementFound)
                            break;
                    }
                }

                if(!replacementFound)
                {
                    log(LOG_LEVEL_VERBOSE, String.format("skipping linked pair(%s) would close existing chain(%d)",
                            newPair.toString(), chain.id()));

                    if(!mSkippedPairs.contains(newPair))
                    {
                        mSkippedPair = true;
                        mSkippedPairs.add(pair);
                    }

                    linkClosesChain = true;
                    continue;
                }
            }

            chain.addLink(newPair, addToStart);
            addedToChain = true;

            LOGGER.debug("index({}) method({}) adding linked pair({} {} len={}) to existing chain({}) {}",
                    mLinkIndex, mLinkReason,
                    newPair.toString(), newPair.assemblyInferredStr(), newPair.length(), chain.id(), addToStart ? "start" : "end");
            break;
        }

        if(!addedToChain)
        {
            if(linkClosesChain)
                return false; // skip this link for now

            SvChain chain = new SvChain(mNextChainId++);
            mChains.add(chain);
            chain.addLink(newPair, true);
            pairToChain[SE_START] = true;
            pairToChain[SE_END] = true;

            LOGGER.debug("index({}) method({}) adding linked pair({} {}) to new chain({})",
                    mLinkIndex, mLinkReason, newPair.toString(), newPair.assemblyInferredStr(), chain.id());
        }

        newPair.setLinkReason(mLinkReason);

        registerNewLink(newPair, pairToChain);
        ++mLinkIndex;

        if(mRunValidation)
            checkHasValidState();

        if(addedToChain)
        {
            // now see if any partial chains can be linked
            reconcileChains();
        }

        return true;
    }

    private SvBreakend findUnlinkedMatchingBreakend(final SvBreakend breakend)
    {
        // get the next available breakend (thereby reducing the replicated instances)
        final List<SvBreakend> breakendList = mUnlinkedBreakendMap.get(breakend.getOrigBreakend());

        if(breakendList == null || breakendList.isEmpty())
            return null;

        return breakendList.get(0);
    }

    private int getUnlinkedBreakendCount(final SvBreakend breakend)
    {
        List<SvBreakend> beList = mUnlinkedBreakendMap.get(breakend);
        return beList != null ? beList.size() : 0;
    }

    private int getMaxUnlinkedPairCount(final SvLinkedPair pair)
    {
        int first = getUnlinkedBreakendCount(pair.getBreakend(true));
        int second = getUnlinkedBreakendCount(pair.getBreakend(false));
        return min(first, second);
    }

    private void addAssemblyLinksToChains()
    {
        if(mAssembledLinks.isEmpty())
            return;

        mLinkReason = "ASMB";

        for(SvLinkedPair pair : mAssembledLinks)
        {
            // add as many times as permitted by the ploidy of the variants
            int pairRepeatCount = !mHasReplication ?
                    1 : min(getSvReplicationCount(pair.first()), getSvReplicationCount(pair.second())) ;

            for(int i = 0; i < pairRepeatCount; ++i)
            {
                addPairToChain(pair);
            }
        }

        if(!mChains.isEmpty())
        {
            LOGGER.debug("created {} partial chains from {} assembly links", mChains.size(), mAssembledLinks.size());
        }
    }

    private void registerNewLink(final SvLinkedPair newPair, boolean[] pairToChain)
    {
        for (int be = SE_START; be <= SE_END; ++be)
        {
            boolean isStart = isStart(be);

            final SvBreakend breakend = newPair.getBreakend(isStart);

            mUnlinkedSVs.remove(breakend.getSV());

            final SvBreakend origBreakend = breakend.getOrigBreakend();

            final List<SvBreakend> breakendList = mUnlinkedBreakendMap.get(origBreakend);

            if(breakendList == null || breakendList.isEmpty())
            {
                LOGGER.error("breakend({}) list already empty", origBreakend.toString());
                mIsValid = false;
                return;
            }

            SvVarData origSV = origBreakend.getSV();
            final SvBreakend otherOrigBreakend = origSV.getBreakend(!breakend.usesStart());

            breakendList.remove(breakend);

            boolean hasUnlinkedBreakend = true;
            if(breakendList.isEmpty())
            {
                mUnlinkedBreakendMap.remove(origBreakend);
                hasUnlinkedBreakend = false;

                // LOGGER.debug("breakend({}) has no more possible links", breakend);

                if(getUnlinkedBreakendCount(otherOrigBreakend) == 0)
                {
                    if(origSV.isFoldback())
                    {
                        // remove if no other instances of this SV remain
                        mFoldbacks.remove(origSV);
                    }
                    else if(mComplexDupCandidates.contains(origSV))
                    {
                        mComplexDupCandidates.remove(origSV);
                    }
                }
            }

            List<SvLinkedPair> possibleLinks = !mSvBreakendPossibleLinks.isEmpty() ? mSvBreakendPossibleLinks.get(origBreakend) : null;

            if(possibleLinks != null)
            {
                if (!hasUnlinkedBreakend)
                {
                    removePossibleLinks(possibleLinks, breakend);
                }
            }

            // check for an opposite pairing between these 2 SVs - need to look into other breakends' lists

            possibleLinks = otherOrigBreakend != null && !mSvBreakendPossibleLinks.isEmpty() ? mSvBreakendPossibleLinks.get(otherOrigBreakend) : null;

            if(possibleLinks != null)
            {
                final SvBreakend otherOrigBreakendAlt = newPair.first() ==  breakend.getSV() ?
                        newPair.second().getOrigSV().getBreakend(!newPair.secondLinkOnStart()) :
                        newPair.first().getOrigSV().getBreakend(!newPair.firstLinkOnStart());

                if(otherOrigBreakendAlt != null)
                {
                    for (SvLinkedPair pair : possibleLinks)
                    {
                        if (pair.hasBreakend(otherOrigBreakend) && pair.hasBreakend(otherOrigBreakendAlt))
                        {
                            possibleLinks.remove(pair);
                            break;
                        }
                    }
                }
            }

            if(mHasReplication)
            {
                // reduce replication counts for breakends which are added to a chain
                if (pairToChain[be])
                {
                    Integer replicationCount = mSvReplicationMap.get(breakend.getOrigSV());

                    if (replicationCount != null)
                    {
                        if (replicationCount <= 2)
                            mSvReplicationMap.remove(breakend.getOrigSV());
                        else
                            mSvReplicationMap.put(breakend.getOrigSV(), replicationCount - 1);
                    }
                }
            }
        }

        // track unique pairs to avoid conflicts (eg end-to-end and start-to-start)
        if(!matchesExistingPair(newPair))
        {
            mUniquePairs.add(newPair);
        }
    }

    private boolean matchesExistingPair(final SvLinkedPair pair)
    {
        for(SvLinkedPair existingPair : mUniquePairs)
        {
            if(pair.matches(existingPair))
                return true;
        }

        return false;
    }

    private void removePossibleLinks(List<SvLinkedPair> possibleLinks, SvBreakend fullyLinkedBreakend)
    {
        if(possibleLinks == null || possibleLinks.isEmpty())
            return;

        final SvVarData linkedSV = fullyLinkedBreakend.getOrigSV();
        final SvBreakend origBreakend = fullyLinkedBreakend.getOrigBreakend();

        int index = 0;
        while (index < possibleLinks.size())
        {
            SvLinkedPair possibleLink = possibleLinks.get(index);
            if (possibleLink.hasBreakend(linkedSV, fullyLinkedBreakend.usesStart()))
            {
                // remove this from consideration
                possibleLinks.remove(index);

                SvBreakend otherBreakend = possibleLink.getBreakend(true) == origBreakend ?
                        possibleLink.getBreakend(false) : possibleLink.getBreakend(true);

                // and remove the pair which was cached in the other breakend's possibles list
                List<SvLinkedPair> otherPossibles = mSvBreakendPossibleLinks.get(otherBreakend);

                if(otherPossibles != null)
                {
                    for (SvLinkedPair otherPair : otherPossibles)
                    {
                        if (otherPair == possibleLink)
                        {
                            otherPossibles.remove(otherPair);

                            if (otherPossibles.isEmpty())
                            {
                                // LOGGER.debug("breakend({}) has no more possible links", otherBreakend);

                                if(!addMorePossibleLinks(otherBreakend, true))
                                {
                                    mSvBreakendPossibleLinks.remove(otherBreakend);
                                }
                            }

                            break;
                        }
                    }
                }
            }
            else
            {
                ++index;
            }
        }

        if(possibleLinks.isEmpty())
        {
            //LOGGER.debug("breakend({}) has no more possible links", origBreakend);

            if(!addMorePossibleLinks(origBreakend, true))
            {
                mSvBreakendPossibleLinks.remove(origBreakend);
            }
        }
    }

    public int getClusterChrBreakendIndex(final SvBreakend breakend)
    {
        if(!mIsClusterSubset)
            return breakend.getClusterChrPosIndex();

        Integer index = mSubsetBreakendClusterIndexMap.get(breakend);
        return index != null ? index : -1;
    }

    private void determinePossibleLinks()
    {
        // form a map of each breakend to its set of all other breakends which can form a valid TI
        // need to exclude breakends which are already assigned to an assembled TI unless replication permits additional instances of it
        // add possible links to a list ordered from shortest to longest length
        // do not chain past a zero cluster allele ploidy
        // identify potential complex DUP candidates along the way
        // for the special case of foldbacks, add every possible link they can make

        /*
        List<SvBreakend> reverseFoldbackBreakends = Lists.newArrayList();

        for(SvVarData foldback : mFoldbacks)
        {
            if(foldback.isChainedFoldback())
            {
                SvBreakend breakend = foldback.getChainedFoldbackBreakend();

                if(breakend.orientation() == 1)
                    reverseFoldbackBreakends.add(breakend);
            }
            else
            {
                if(foldback.orientation(true) == 1)
                {
                    reverseFoldbackBreakends.add(foldback.getBreakend(true));
                    reverseFoldbackBreakends.add(foldback.getBreakend(false));
                }
            }
        }
        */

        for (final Map.Entry<String, List<SvBreakend>> entry : mChrBreakendMap.entrySet())
        {
            final String chromosome = entry.getKey();
            final List<SvBreakend> breakendList = entry.getValue();
            final double[][] allelePloidies = mClusterPloidyLimits.getChrAllelePloidies().get(chromosome);

            for (int i = 0; i < breakendList.size() -1; ++i)
            {
                final SvBreakend lowerBreakend = breakendList.get(i);

                boolean matchedPloidy = false;

                if(lowerBreakend.orientation() != -1)
                    continue;

                if(alreadyLinkedBreakend(lowerBreakend))
                    continue;

                List<SvLinkedPair> lowerPairs = mSvBreakendPossibleLinks.get(lowerBreakend);

                if(lowerPairs == null)
                {
                    lowerPairs = Lists.newArrayList();
                    mSvBreakendPossibleLinks.put(lowerBreakend, lowerPairs);
                }

                final SvVarData lowerSV = lowerBreakend.getSV();

                boolean lowerValidAP = mUseAllelePloidies && mClusterPloidyLimits.hasValidAllelePloidyData(
                        getClusterChrBreakendIndex(lowerBreakend), allelePloidies);

                boolean lowerIsFoldback = lowerSV.isFoldback() && (!lowerSV.isChainedFoldback() || lowerSV.getChainedFoldbackBreakend() == lowerBreakend);

                int skippedNonAssembledIndex = -1; // the first index of a non-assembled breakend after the current one

                for (int j = i+1; j < breakendList.size(); ++j)
                {
                    final SvBreakend upperBreakend = breakendList.get(j);

                    if(skippedNonAssembledIndex == -1)
                    {
                        if(!upperBreakend.isAssembledLink())
                        {
                            // invalidate the possibility of these 2 breakends satisfying the complex DUP scenario
                            skippedNonAssembledIndex = j;
                        }
                    }

                    if(upperBreakend.orientation() != 1)
                        continue;

                    if(upperBreakend.getSV() == lowerBreakend.getSV())
                        continue;

                    if(alreadyLinkedBreakend(upperBreakend))
                        continue;

                    long distance = upperBreakend.position() - lowerBreakend.position();
                    int minTiLength = getMinTemplatedInsertionLength(lowerBreakend, upperBreakend);

                    if(distance < minTiLength)
                        continue;

                    // record the possible link
                    final SvVarData upperSV = upperBreakend.getSV();

                    SvLinkedPair newPair = new SvLinkedPair(lowerSV, upperSV, LINK_TYPE_TI,
                            lowerBreakend.usesStart(), upperBreakend.usesStart());

                    // if(getSvReplicationCount(lowerSV) == getSvReplicationCount(upperSV))
                        // matchedPloidy = true;

                    // make note of any pairs formed from adjacent facing breakends
                    if(j == i + 1 && getSvReplicationCount(lowerSV) == getSvReplicationCount(upperSV))
                    {
                        mAdjacentMatchingPairs.add(newPair);
                    }

                    lowerPairs.add(newPair);

                    List<SvLinkedPair> upperPairs = mSvBreakendPossibleLinks.get(upperBreakend);

                    if(upperPairs == null)
                    {
                        upperPairs = Lists.newArrayList();
                        mSvBreakendPossibleLinks.put(upperBreakend, upperPairs);

                        // create an entry at the upper breakend's start point to indicate it hasn't begun its search
                        // mBreakendLastLinkIndexMap.put(upperBreakend, getClusterChrBreakendIndex(upperBreakend));
                    }

                    upperPairs.add(0, newPair); // add to front since always nearer than the one prior

                    if(skippedNonAssembledIndex == -1 || skippedNonAssembledIndex == j)
                    {
                        // make note of any breakends which run into a high-ploidy SV at their first opposing breakend
                        if (!lowerBreakend.getSV().isFoldback())
                        {
                            checkIsComplexDupSV(lowerBreakend, upperBreakend);
                        }

                        if (!upperBreakend.getSV().isFoldback())
                        {
                            checkIsComplexDupSV(upperBreakend, lowerBreakend);
                        }
                    }

                    if(lowerValidAP && mClusterPloidyLimits.hasValidAllelePloidyData(
                            getClusterChrBreakendIndex(upperBreakend), allelePloidies))
                    {
                        double clusterAP = allelePloidies[getClusterChrBreakendIndex(upperBreakend)][CLUSTER_AP];

                        if(clusterAP < CLUSTER_ALLELE_PLOIDY_MIN)
                        {
                            // this lower breakend cannot match with anything further upstream
                            log(LOG_LEVEL_VERBOSE, String.format("breakends lower(%d: %s) limited at upper(%d: %s) with clusterAP(%.2f)",
                                    i, lowerBreakend.toString(), j, upperBreakend.toString(), clusterAP));

                            break;
                        }
                    }

                    /* currently disabled
                    if(matchedPloidy && exceedsMaxPossibleLinks(lowerPairs.size()) && !lowerIsFoldback)
                    {
                        // more possible links could be craeted, but pause adding any more for now
                        mBreakendLastLinkIndexMap.put(lowerBreakend, j);
                        break;
                    }
                    */
                }
            }
        }

        /*
        for(SvBreakend breakend : reverseFoldbackBreakends)
        {
            mBreakendLastLinkIndexMap.put(breakend, getClusterChrBreakendIndex(breakend));
            addMorePossibleLinks(breakend, false);
        }
        */
    }

    private void checkIsComplexDupSV(SvBreakend lowerPloidyBreakend, SvBreakend higherPloidyBreakend)
    {
        SvVarData var = lowerPloidyBreakend.getSV();

        if(var.isNullBreakend() || var.type() == DEL)
            return;

        if(mComplexDupCandidates.contains(var))
            return;

        if(var.ploidyMin() * 2 > higherPloidyBreakend.getSV().ploidyMax())
            return;

        boolean lessThanMax = var.ploidyMax() < higherPloidyBreakend.getSV().ploidyMin();

        // check whether the other breakend satisfies the same ploidy comparison criteria
        SvBreakend otherBreakend = var.getBreakend(!lowerPloidyBreakend.usesStart());

        final List<SvBreakend> breakendList = mChrBreakendMap.get(otherBreakend.chromosome());

        boolean traverseUp = otherBreakend.orientation() == -1;
        int index = getClusterChrBreakendIndex(otherBreakend);

        while(true)
        {
            index += traverseUp ? 1 : -1;

            if(index < 0 || index >= breakendList.size())
                break;

            final SvBreakend breakend = breakendList.get(index);

            if(breakend == lowerPloidyBreakend)
                break;

            if (breakend.isAssembledLink())
            {
                index += traverseUp ? 1 : -1;
                continue;
            }

            if (breakend.orientation() == otherBreakend.orientation())
                break;

            SvVarData otherSV = breakend.getSV();

            if(var.ploidyMin() * 2 <= otherSV.ploidyMax())
            {
                if(lessThanMax || var.ploidyMax() < otherSV.ploidyMin())
                {
                    LOGGER.debug("identified complex dup({} {}) vs SVs({} & {})",
                            var.posId(), var.type(), higherPloidyBreakend.getSV().id(), otherSV.id());
                    mComplexDupCandidates.add(var);
                }
            }

            break;
        }
    }

    private boolean alreadyLinkedBreakend(final SvBreakend breakend)
    {
        // assembled links have already been added to chains prior to determining remaining possible links
        // so these need to be excluded unless their replication count allows them to be used again
        return breakend.isAssembledLink() && getUnlinkedBreakendCount(breakend) == 0;
    }

    private int getSvReplicationCount(final SvVarData var)
    {
        Integer repCount = mSvOriginalReplicationMap.get(var);
        return repCount != null ? repCount.intValue() : 1;
    }

    private void setSvReplicationCounts()
    {
        if(!mHasReplication)
            return;

        for(final SvVarData var : mUniqueSVs)
        {
            if(var.getReplicatedCount() > 0)
            {
                mSvReplicationMap.put(var, var.getReplicatedCount());
                mSvOriginalReplicationMap.put(var, var.getReplicatedCount());
            }
            else
            {
                mSvOriginalReplicationMap.put(var, 1);
            }
        }
    }

    private void setUnlinkedBreakends()
    {
        // make a cache of all unchained breakends in those of replicated SVs
        for(final SvVarData var : mReplicatedSVs)
        {
            for (int be = SE_START; be <= SE_END; ++be)
            {
                boolean isStart = isStart(be);

                if (var.isNullBreakend() && !isStart)
                    continue;

                final SvBreakend breakend = var.getBreakend(isStart);
                final SvBreakend origBreakend = breakend.getOrigBreakend();

                List<SvBreakend> breakends = mUnlinkedBreakendMap.get(origBreakend);

                if(breakends == null)
                {
                    breakends = Lists.newArrayList();
                    mUnlinkedBreakendMap.put(origBreakend, breakends);
                }

                breakends.add(breakend);
            }
        }

        mUnlinkedSVs.addAll(mReplicatedSVs);
    }

    private List<SvVarData> getMaxReplicationSVs()
    {
        List<SvVarData> maxRepIds = Lists.newArrayList();
        int maxRepCount = 1;

        for(Map.Entry<SvVarData,Integer> entry : mSvReplicationMap.entrySet())
        {
            int repCount = entry.getValue();

            if(repCount <= 1)
                continue;

            if(repCount > maxRepCount)
            {
                maxRepCount = repCount;
                maxRepIds.clear();
                maxRepIds.add(entry.getKey());
            }
            else if(repCount == maxRepCount)
            {
                if(!maxRepIds.contains(entry.getKey()))
                    maxRepIds.add(entry.getKey());
            }
        }

        return maxRepIds;
    }

    private void reconcileChains()
    {
        int index1 = 0;
        while(index1 < mChains.size())
        {
            SvChain chain1 = mChains.get(index1);

            boolean chainsMerged = false;

            for (int index2 = index1 + 1; index2 < mChains.size(); ++index2)
            {
                SvChain chain2 = mChains.get(index2);

                for (int be1 = SE_START; be1 <= SE_END; ++be1)
                {
                    boolean c1Start = isStart(be1);

                    for (int be2 = SE_START; be2 <= SE_END; ++be2)
                    {
                        boolean c2Start = isStart(be2);

                        if (chain1.canAddLinkedPair(chain2.getLinkedPair(c2Start), c1Start, false))
                        {
                            LOGGER.debug("merging chain({} links={}) with chain({} links={})",
                                    chain1.id(), chain1.getLinkCount(), chain2.id(), chain2.getLinkCount());

                            if(c2Start)
                            {
                                // merge chains and remove the latter
                                for (SvLinkedPair linkedPair : chain2.getLinkedPairs())
                                {
                                    chain1.addLink(linkedPair, c1Start);
                                }
                            }
                            else
                            {
                                // add in reverse
                                for (int index = chain2.getLinkedPairs().size() - 1; index >= 0; --index)
                                {
                                    SvLinkedPair linkedPair = chain2.getLinkedPairs().get(index);
                                    chain1.addLink(linkedPair, c1Start);
                                }
                            }

                            mChains.remove(index2);

                            chainsMerged = true;
                            break;
                        }

                    }

                    if (chainsMerged)
                        break;
                }

                if (chainsMerged)
                    break;
            }

            if (!chainsMerged)
            {
                ++index1;
            }
        }
    }

    private static int LOG_LEVEL_ERROR = 0;
    private static int LOG_LEVEL_INFO = 1;
    private static int LOG_LEVEL_DEBUG = 2;
    private static int LOG_LEVEL_VERBOSE = 3;

    private void log(int level, final String message)
    {
        if(level >= LOG_LEVEL_VERBOSE && !mLogVerbose)
            return;

        if(level >= LOG_LEVEL_DEBUG && !LOGGER.isDebugEnabled())
            return;

        LOGGER.debug(message);
    }

    private void checkProgress()
    {
        if(!LOGGER.isDebugEnabled())
            return;

        if(!mHasReplication || mUniqueSVs.size() < 100)
            return;

        if((mLinkIndex % 100) == 0)
        {
            LOGGER.debug("cluster({}) chaining progress: SVs({}) partialChains({}) unlinked(SVs={} breakends={}) replicatedSVs({})",
                    mClusterId, mUniqueSVs.size(), mChains.size(), mUnlinkedSVs.size(),
                    mUnlinkedBreakendMap.size(), mSvReplicationMap.size());
        }
    }

    private boolean checkHasValidState()
    {
        // first check that the remaining possible links are supported by unlinked breakends
        for(Map.Entry<SvBreakend,List<SvLinkedPair>> entry : mSvBreakendPossibleLinks.entrySet())
        {
            SvBreakend breakend = entry.getKey();

            List<SvBreakend> breakendList = mUnlinkedBreakendMap.get(breakend);

            if(breakendList == null)
            {
                LOGGER.error("cluster({}) runIndex({}): breakend({}) has {} possible pairs but no available breakends",
                        mClusterId, mLinkIndex, breakend.toString(), entry.getValue().size());

                mIsValid = false;
            }
        }

        return mIsValid;
    }

    public void setMaxPossibleLinks(int maxLinks)
    {
        LOGGER.warn("incremental link finding disabled");

        /*
        if(maxLinks == 0)
            mMaxPossibleLinks = 0;
        else
            mMaxPossibleLinks = max(maxLinks, 2);
        */
    }

    private boolean isOppositeMatchVsExisting(final SvLinkedPair pair)
    {
        for(SvLinkedPair existingPair : mUniquePairs)
        {
            if(pair.oppositeMatch(existingPair))
                return true;
        }

        return false;
    }

    private boolean exceedsMaxPossibleLinks(int linkCount)
    {
        return mMaxPossibleLinks > 0 && linkCount >= mMaxPossibleLinks;
    }

    private boolean addMorePossibleLinks(SvBreakend breakend, boolean applyMax)
    {
        return false;

        /*
        Integer lastIndex = mBreakendLastLinkIndexMap.get(breakend);

        if(lastIndex == null || lastIndex < 0)
            return false;

        if(getUnlinkedBreakendCount(breakend) == 0)
        {
            mBreakendLastLinkIndexMap.remove(breakend);
            return false;
        }

        // begin from immediately after the last added index and try to add another X possible links
        final List<SvBreakend> breakendList = mChrBreakendMap.get(breakend.chromosome());

        final double[][] allelePloidies = mChrAllelePloidies.get(breakend.chromosome());

        boolean hasValidAP = mUseAllelePloidies && hasValidAllelePloidyData(breakend, allelePloidies);

        List<SvLinkedPair> possiblePairs = mSvBreakendPossibleLinks.get(breakend);

        if(possiblePairs == null)
            return false;

        boolean traverseUp = breakend.orientation() == -1;
        int index = lastIndex;

        boolean matchedPloidy = false;
        int linksAdded = 0;
        boolean lastIndexValid = true;
        while (true)
        {
            index += traverseUp ? 1 : -1;

            if(index < 0 || index >= breakendList.size())
            {
                lastIndexValid = false;
                break;
            }

            final SvBreakend otherBreakend = breakendList.get(index);

            if(otherBreakend.orientation() == breakend.orientation())
                continue;

            if(otherBreakend.getSV() == breakend.getSV())
                continue;

            if(getUnlinkedBreakendCount(otherBreakend) == 0)
                continue;

            long distance = abs(otherBreakend.position() - breakend.position());
            int minTiLength = getMinTemplatedInsertionLength(breakend, otherBreakend);

            if(distance < minTiLength)
                continue;

            List<SvLinkedPair> otherPairs = mSvBreakendPossibleLinks.get(otherBreakend);

            if(otherPairs == null)
                continue;

            // record the possible link
            SvBreakend lowerBreakend = breakend.orientation() == -1 ? breakend : otherBreakend;
            SvBreakend upperBreakend = breakend.orientation() == 1 ? breakend : otherBreakend;
            final SvVarData lowerSV = lowerBreakend.getSV();
            final SvVarData upperSV = upperBreakend.getSV();

            SvLinkedPair newPair = new SvLinkedPair(lowerSV, upperSV, LINK_TYPE_TI,
                    lowerBreakend.usesStart(), upperBreakend.usesStart());

            // check link hasn't already been added (which can happen if added from the other breakend)
            boolean skipPair = false;

            for(SvLinkedPair existingPair : possiblePairs)
            {
                if(existingPair.matches(newPair))
                {
                    skipPair = true;
                    break;
                }
            }

            if(skipPair)
                continue;

            for(SvLinkedPair existingPair : otherPairs)
            {
                if(existingPair.matches(newPair))
                {
                    skipPair = true;
                    break;
                }
            }

            if(skipPair)
                continue;

            // check for a clash against existing pairs
            if(isOppositeMatchVsExisting(newPair))
                continue;

            ++linksAdded;
            possiblePairs.add(newPair);
            otherPairs.add(newPair);

            if(!matchedPloidy)
            {
                matchedPloidy = getSvReplicationCount(lowerSV) == getSvReplicationCount(upperSV);
            }

            if(hasValidAP && hasValidAllelePloidyData(otherBreakend, allelePloidies))
            {
                double clusterAP = allelePloidies[getClusterChrBreakendIndex(otherBreakend)][CLUSTER_AP];

                if(clusterAP < CLUSTER_ALLELE_PLOIDY_MIN)
                {
                    // this lower breakend cannot match with anything further upstream
                    log(LOG_LEVEL_VERBOSE, String.format("breakend(%d: %s) limited by other(%d: %s) with clusterAP(%.2f)",
                            getClusterChrBreakendIndex(breakend), breakend.toString(), index, otherBreakend.toString(), clusterAP));

                    lastIndexValid = false;
                    break;
                }
            }

            if(applyMax && matchedPloidy && exceedsMaxPossibleLinks(possiblePairs.size()))
            {
                break;
            }
        }

        if(lastIndexValid)
        {
            // make note of the last location tested for adding a new possible link
            mBreakendLastLinkIndexMap.put(breakend, index);
        }
        else
        {
            mBreakendLastLinkIndexMap.remove(breakend);
        }

        return linksAdded > 0;
        */
    }

    private void cullPossibleLinks()
    {
        if(mMaxPossibleLinks == 0)
            return;

        int culledPairs = 0;
        for(Map.Entry<SvBreakend,List<SvLinkedPair>> entry : mSvBreakendPossibleLinks.entrySet())
        {
            SvBreakend breakend = entry.getKey();

            List<SvLinkedPair> possiblePairs = mSvBreakendPossibleLinks.get(breakend);

            if(exceedsMaxPossibleLinks(possiblePairs.size()))
                continue;

            int index = possiblePairs.size() - 1;
            while(exceedsMaxPossibleLinks(index))
            {
                SvLinkedPair pair = possiblePairs.get(index);
                SvBreakend otherBreakend = pair.getOtherBreakend(breakend);

                // only remove if the other breakend also has an excess of possible pairs and it's not in the first X entries
                List<SvLinkedPair> otherPossiblePairs = mSvBreakendPossibleLinks.get(otherBreakend);

                boolean canRemoveOther = true;
                for(int index2 = 0; index2 < mMaxPossibleLinks; ++index2)
                {
                    if(otherPossiblePairs.get(index2) == pair)
                    {
                        canRemoveOther = false;
                        break;
                    }
                }

                if(canRemoveOther)
                {
                    otherPossiblePairs.remove(pair);
                    possiblePairs.remove(possiblePairs.size() - 1);
                    ++culledPairs;
                }

                --index;
            }
        }

        LOGGER.debug("cluster({}) culled {} possible pairs", mClusterId, culledPairs);
    }

}
