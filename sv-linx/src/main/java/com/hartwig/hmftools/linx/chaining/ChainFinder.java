package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatPloidy;
import static com.hartwig.hmftools.linx.chaining.ChainDiagnostics.LOG_TYPE_INFO;
import static com.hartwig.hmftools.linx.chaining.ChainDiagnostics.LOG_TYPE_VERBOSE;
import static com.hartwig.hmftools.linx.chaining.ChainDiagnostics.LOG_TYPE_WARN;
import static com.hartwig.hmftools.linx.chaining.ChainPloidyLimits.CLUSTER_ALLELE_PLOIDY_MIN;
import static com.hartwig.hmftools.linx.chaining.ChainPloidyLimits.CLUSTER_AP;
import static com.hartwig.hmftools.linx.chaining.ChainPloidyLimits.calcPloidyUncertainty;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.areLinkedSection;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.getMinTemplatedInsertionLength;
import static com.hartwig.hmftools.linx.chaining.SvChain.checkIsValid;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LINK_TYPE_TI;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.cn.PloidyCalcData;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

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
    public String mSampleId;
    private boolean mHasReplication;

    // input state
    private final List<SvVarData> mSvList;
    private final List<SvVarData> mFoldbacks;
    private final List<SvVarData> mDoubleMinuteSVs;
    private final List<SvLinkedPair> mAssembledLinks;
    private Map<String,List<SvBreakend>> mChrBreakendMap;

    // links a breakend to its position in a chromosome's breakend list - only applicable if a subset of SVs are being chained
    private boolean mIsClusterSubset;
    private final Map<SvBreakend,Integer> mSubsetBreakendClusterIndexMap;

    // chaining state
    private final List<SvLinkedPair> mSkippedPairs; // pairs which would close a chain, temporarily ignored
    private final List<SvVarData> mComplexDupCandidates; // identified SVs which duplication another SV
    private final List<SvLinkedPair> mUniquePairs; // cache of unique pairs added through c
    private final List<SvLinkedPair> mAdjacentMatchingPairs;
    private boolean mPairSkipped; // keep track of any excluded pair or SV without exiting the chaining routine

    private final List<SvChain> mChains;
    private final List<SvChain> mUniqueChains;
    private int mNextChainId;

    // a cache of cluster ploidy boundaries which links cannot cross
    private final ChainPloidyLimits mClusterPloidyLimits;

    // determined up-front - the set of all possible links from a specific breakend to other breakends
    private final Map<SvBreakend, List<SvLinkedPair>> mSvBreakendPossibleLinks;

    //
    private final Map<SvVarData, SvChainState> mSvConnectionsMap;
    private final List<SvChainState> mSvCompletedConnections;
    private List<SvVarData> mReplicatedSVs;
    private List<SvBreakend> mReplicatedBreakends;

    // temporary support for old chain finder
    private ChainFinderOld mOldFinder;
    private boolean mUseOld;

    public static final int CHAIN_METHOD_OLD = 0;
    public static final int CHAIN_METHOD_NEW = 1;

    private int mLinkIndex; // incrementing value for each link added to any chain
    private String mLinkReason;
    private boolean mIsValid;
    private boolean mLogVerbose;
    private boolean mRunValidation;
    private boolean mUseAllelePloidies;

    private static final String LR_METHOD_ASMB = "ASMB";
    private static final String LR_METHOD_ONLY = "ONLY";
    private static final String LR_METHOD_FOLDBACK = "FOLDBACK";
    private static final String LR_METHOD_ADJAC = "ADJAC";
    private static final String LR_METHOD_PL_MATCH = "PL_MATCH";
    private static final String LR_METHOD_PL_MAX = "PL_MAX";
    private static final String LR_METHOD_SHORT = "SHORT";
    private static final String LR_METHOD_CMP_DUP = "CMP_DUP";
    private static final String LR_METHOD_DM_DUP = "DM_DUP";
    private static final String LR_METHOD_DM_CLOSE = "DM_CLOSE";

    // self-analysis only
    private final ChainDiagnostics mDiagnostics;

    private static final Logger LOGGER = LogManager.getLogger(ChainFinder.class);

    public ChainFinder()
    {
        mSvList = Lists.newArrayList();
        mFoldbacks = Lists.newArrayList();
        mDoubleMinuteSVs = Lists.newArrayList();
        mIsClusterSubset = false;
        mAssembledLinks = Lists.newArrayList();
        mChrBreakendMap = null;

        mClusterPloidyLimits = new ChainPloidyLimits();

        mAdjacentMatchingPairs = Lists.newArrayList();
        mSvCompletedConnections = Lists.newArrayList();
        mSubsetBreakendClusterIndexMap = Maps.newHashMap();
        mSvConnectionsMap = Maps.newHashMap();
        mComplexDupCandidates = Lists.newArrayList();
        mChains = Lists.newArrayList();
        mUniqueChains = Lists.newArrayList();
        mSkippedPairs = Lists.newArrayList();
        mSvBreakendPossibleLinks = Maps.newHashMap();
        mUniquePairs = Lists.newArrayList();
        mReplicatedSVs = Lists.newArrayList();
        mReplicatedBreakends = Lists.newArrayList();

        mHasReplication = false;
        mLogVerbose = false;
        mRunValidation = false;
        mIsValid = true;
        mPairSkipped = false;
        mNextChainId = 0;
        mLinkIndex = 0;
        mLinkReason = "";
        mSampleId= "";
        mUseAllelePloidies = false;

        mDiagnostics = new ChainDiagnostics(
                mSvConnectionsMap, mSvCompletedConnections, mChains, mUniqueChains,
                mSvBreakendPossibleLinks, mDoubleMinuteSVs, mUniquePairs);

        mOldFinder = new ChainFinderOld();
        mUseOld = false;
    }

    public void setUseOldMethod(boolean toggle)
    {
        LOGGER.info("using {} chain-finder", toggle ? "old" : "new");
        mUseOld = toggle;
    }

    public void clear()
    {
        mClusterId = -1;
        mSvList.clear();
        mSvConnectionsMap.clear();
        mSvCompletedConnections.clear();
        mFoldbacks.clear();
        mDoubleMinuteSVs.clear();
        mIsClusterSubset = false;
        mAssembledLinks.clear();
        mChrBreakendMap = null;

        mAdjacentMatchingPairs.clear();
        mSubsetBreakendClusterIndexMap.clear();
        mComplexDupCandidates.clear();
        mChains.clear();
        mUniqueChains.clear();
        mSkippedPairs.clear();
        mSvBreakendPossibleLinks.clear();
        mUniquePairs.clear();
        mReplicatedSVs.clear();
        mReplicatedBreakends.clear();

        mNextChainId = 0;
        mLinkIndex = 0;
        mIsValid = true;
        mPairSkipped = false;

        mDiagnostics.clear();
    }

    public void setSampleId(final String sampleId)
    {
        mSampleId = sampleId;
        mDiagnostics.setSampleId(sampleId);
    }

    public void initialise(SvCluster cluster)
    {
        // attempt to chain all the SVs in a cluster

        // critical that all state is cleared before the next run
        clear();

        // isSpecificCluster(cluster);

        mClusterId = cluster.id();
        mSvList.addAll(cluster.getSVs());
        mFoldbacks.addAll(cluster.getFoldbacks());
        mDoubleMinuteSVs.addAll(cluster.getDoubleMinuteSVs());
        mAssembledLinks.addAll(cluster.getAssemblyLinkedPairs());
        mChrBreakendMap = cluster.getChrBreakendMap();
        mHasReplication = cluster.requiresReplication();
        mIsClusterSubset = false;

        if(mUseOld)
            mOldFinder.initialise(cluster);
    }

    public void initialise(SvCluster cluster, final List<SvVarData> svList)
    {
        // chain a specific subset of a cluster's SVs
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
            if(!mHasReplication && var.isReplicatedSv())
            {
                mHasReplication = true;
                continue;
            }

            mSvList.add(var);

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

        mDoubleMinuteSVs.addAll(cluster.getDoubleMinuteSVs());

        if(mUseOld)
            mOldFinder.initialise(cluster, svList);
    }

    public void setLogVerbose(boolean toggle)
    {
        mLogVerbose = toggle;
        setRunValidation(toggle);

        if(mUseOld)
            mOldFinder.setLogVerbose(toggle);
    }

    public void setRunValidation(boolean toggle) { mRunValidation = toggle; }
    public void setUseAllelePloidies(boolean toggle) { mUseAllelePloidies = toggle; }

    public final List<SvChain> getUniqueChains()
    {
        return mUseOld ? mOldFinder.getUniqueChains() : mUniqueChains;
    }
    public double getValidAllelePloidySegmentPerc() { return mClusterPloidyLimits.getValidAllelePloidySegmentPerc(); }
    public final ChainDiagnostics getDiagnostics() { return mDiagnostics; }

    public void formChains(boolean assembledLinksOnly)
    {
        if(mUseOld)
        {
            mOldFinder.formChains(assembledLinksOnly);
            return;
        }

        if(mSvList.size() < 2)
            return;

        if (mSvList.size() >= 4)
        {
            LOGGER.debug("cluster({}) starting chaining with assemblyLinks({}) svCount({})",
                    mClusterId, mAssembledLinks.size(), mSvList.size());
        }

        mClusterPloidyLimits.initialise(mClusterId, mChrBreakendMap);

        buildChains(assembledLinksOnly);

        checkChains();
        removeIdenticalChains();

        // TODO add any replicated SVs to the cluster
        // mReplicatedSVs

        mDiagnostics.chainingComplete();

        if(!mIsValid)
        {
            LOGGER.warn("cluster({}) chain finding failed", mClusterId);
            return;
        }
    }

    private void checkChains()
    {
        for(final SvChain chain : mChains)
        {
            if(!checkIsValid(chain))
            {
                LOGGER.error("cluster({}) has invalid chain({})");
                chain.logLinks();
                mIsValid = false;
            }
        }
    }

    private void removeIdenticalChains()
    {
        if(!mHasReplication)
        {
            mUniqueChains.addAll(mChains);
            return;
        }

        for(final SvChain newChain : mChains)
        {
            boolean matched = false;

            for(final SvChain chain : mUniqueChains)
            {
                if (chain.identicalChain(newChain, true))
                {
                    LOGGER.debug("cluster({}) skipping duplicate chain({}) vs origChain({})",
                            mClusterId, newChain.id(), chain.id());

                    // record repeated links
                    for(SvLinkedPair pair : chain.getLinkedPairs())
                    {
                        if(newChain.getLinkedPairs().stream().anyMatch(x -> x.matches(pair)))
                        {
                            pair.setRepeatCount(pair.repeatCount()+1);
                        }
                    }

                    matched = true;
                    break;
                }
            }

            if(!matched)
            {
                mUniqueChains.add(newChain);
            }
        }
    }

    public void addChains(SvCluster cluster)
    {
        if(mUseOld)
        {
            mOldFinder.addChains(cluster);
            return;
        }

        // add these chains to the cluster, but skip any which are identical to existing ones,
        // which can happen for clusters with replicated SVs
        mUniqueChains.stream().forEach(chain -> checkAddNewChain(chain, cluster));

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
        if(!mHasReplication)
        {
            cluster.addChain(newChain, false, true);
            return;
        }

        // any identical chains (including precise subsets of longer chains) will have their replicated SVs entirely removed
        if(!mUniqueChains.contains(newChain))
        {
            /* the replicated SVs from outside the chain finder will no longer match
            boolean allReplicatedSVs = newChain.getSvCount(false) == 0;

            // remove these replicated SVs as well as the replicated chain
            if(allReplicatedSVs)
            {
                newChain.getSvList().stream().forEach(x -> cluster.removeReplicatedSv(x));
            }
            */

            return;
        }

        cluster.addChain(newChain, false, true);
    }

    private void buildChains(boolean assembledLinksOnly)
    {
        populateSvPloidyMap();

        mDiagnostics.initialise(mClusterId, mHasReplication);

        // first make chains out of any assembly links
        addAssemblyLinksToChains();

        if(assembledLinksOnly)
            return;

        if(mUseAllelePloidies && mHasReplication)
            mClusterPloidyLimits.determineBreakendPloidies();

        determinePossibleLinks();

        mDiagnostics.setPriorityData(mComplexDupCandidates, mFoldbacks);

        int iterationsWithoutNewLinks = 0; // protection against loops

        while (true)
        {
            mPairSkipped = false;
            int lastAddedIndex = mLinkIndex;

            List<ProposedLinks> proposedLinks = findProposedLinks();

            if(proposedLinks.isEmpty())
            {
                if(!mPairSkipped)
                    break;
            }
            else
            {
                processProposedLinks(proposedLinks);

                if(!mIsValid)
                    return;
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

            mDiagnostics.checkProgress(mLinkIndex);
        }

        checkDoubleMinuteChains();
    }

    private List<ProposedLinks> findProposedLinks()
    {
        // proceed through the link-finding methods in priority
        List<ProposedLinks> proposedLinks = findSingleOptionPairs();

        mLinkReason = "";

        if(proposedLinks != null && !proposedLinks.isEmpty())
        {
            mLinkReason = LR_METHOD_ONLY;
            return proposedLinks;
        }

        if (mHasReplication)
        {
            proposedLinks = findFoldbackPairs();

            if (proposedLinks != null)
            {
                mLinkReason = LR_METHOD_FOLDBACK;
                return proposedLinks;
            }
        }

        proposedLinks = findAdjacentMatchingPairs();
        if (proposedLinks != null)
        {
            mLinkReason = LR_METHOD_ADJAC;
            return proposedLinks;
        }

        if(mHasReplication)
        {
            proposedLinks = findPloidyMatchPairs();

            if (proposedLinks != null)
            {
                mLinkReason = LR_METHOD_PL_MATCH;
                return proposedLinks;
            }

            proposedLinks = findComplexDupPairs();

            if (proposedLinks != null)
            {
                mLinkReason = LR_METHOD_CMP_DUP;
                return proposedLinks;
            }

            // proposedLinks = findMaxReplicationPairs();

            if (proposedLinks != null)
            {
                mLinkReason = LR_METHOD_PL_MAX;
                return proposedLinks;
            }
        }

        // try all remaining

        //List<SvBreakend> breakendList = mUnlinkedBreakendMap.keySet().stream().collect(Collectors.toList());

        List<SvBreakend> breakendList = Lists.newArrayList();

        for(SvChainState svConn : mSvConnectionsMap.values())
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(svConn.SV.isNullBreakend() && se == SE_END)
                    continue;

                if (svConn.unlinked(se) > 0)
                    breakendList.add(svConn.SV.getBreakend(isStart(se)));
            }
        }

        // proposedLinks = findFewestOptionPairs(breakendList, false);

        if(proposedLinks != null)
        {
            mLinkReason = LR_METHOD_SHORT;
            return proposedLinks;
        }

        return proposedLinks;
    }

    private List<ProposedLinks> findPloidyMatchPairs()
    {
        ProposedLinks bestProposedLink = null;

        // find the highest pair of SVs with matching ploidy
        double maxPloidyMatch = 0;

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

                if(breakendPloidy >= maxPloidyMatch)
                {
                    SvBreakend breakend = var.getBreakend(isStart);
                    List<SvLinkedPair> svLinks = mSvBreakendPossibleLinks.get(breakend);

                    if(svLinks == null)
                        continue;

                    for(SvLinkedPair pair : svLinks)
                    {
                        if(bestProposedLink != null && bestProposedLink.Links.contains(pair))
                            continue;

                        if(mSkippedPairs.contains(pair))
                            continue;

                        SvBreakend otherBreakend = pair.getOtherBreakend(breakend);
                        SvChainState otherConn = mSvConnectionsMap.get(otherBreakend.getSV());

                        if(otherConn == null)
                            continue;

                        double otherBreakendPloidy = otherConn.unlinked(otherBreakend.usesStart());

                        if(!copyNumbersEqual(otherBreakendPloidy, breakendPloidy))
                            continue;

                        SvVarData otherVar = otherBreakend.getSV();

                        log(LOG_TYPE_VERBOSE, String.format("pair(%s) with matching ploidy(%s & %s)",
                                pair.toString(), formatPloidy(breakendPloidy), formatPloidy(otherBreakendPloidy)));

                        double avgPloidy = (breakendPloidy + otherBreakendPloidy) * 0.5;

                        if(bestProposedLink == null || avgPloidy > bestProposedLink.Ploidy)
                        {
                            // take if higher
                            bestProposedLink = new ProposedLinks(pair, avgPloidy);
                        }
                        else
                        {
                            continue;
                        }

                        if(var.getImpliedPloidy() != otherVar.getImpliedPloidy())
                        {
                            log(LOG_TYPE_WARN, String.format("ploidy-match SVs(%s & %s) have diff ploidies(%d & %d)",
                                    var.id(), otherVar.id(), var.getImpliedPloidy(), otherVar.getImpliedPloidy()));
                        }
                    }
                }
            }
        }

        return bestProposedLink != null ? Lists.newArrayList(bestProposedLink) : Lists.newArrayList();
    }

    private List<ProposedLinks> findComplexDupPairs()
    {
        // both ends of a foldback or complex DUP connect to one end of another SV with ploidy >= 2x
        if(mComplexDupCandidates.isEmpty())
            return Lists.newArrayList();

        List<SvLinkedPair> bestCompDupLinks = Lists.newArrayList();
        double maxCompDupPloidy = 0;
        SvChain matchedChain = null;

        // the complex DUP SVs need to connect to both ends of either a single SV or a set of chained SVs with twice the ploidy
        // for a complex DUP of the form D - A - D, where A has double ploidy of D, and both ends of D connect to both ends of A
        for(SvVarData compDup : mComplexDupCandidates)
        {
            int compDupPloidy = compDup.getImpliedPloidy();

            SvBreakend compDupBeStart = compDup.getBreakend(true);
            SvBreakend compDupBeEnd = compDup.getBreakend(false);

            List<SvLinkedPair> pairsOnCdStart = mSvBreakendPossibleLinks.get(compDupBeStart);
            List<SvLinkedPair> pairsOnCdEnd = mSvBreakendPossibleLinks.get(compDupBeEnd);

            if (pairsOnCdStart == null || pairsOnCdEnd == null)
                continue;

            List<SvLinkedPair[]> candidatePairs = Lists.newArrayList();

            // search existing chains for open chain ends match the set of possibles for the complex DUP and with twice the ploidy
            for(SvChain chain : mChains)
            {
                SvBreakend chainBeStart = chain.getOpenBreakend(true);
                SvBreakend chainBeEnd = chain.getOpenBreakend(false);

                if(chainBeStart == null || chainBeEnd == null)
                    continue;

                chainBeStart = chainBeStart.getOrigBreakend();
                chainBeEnd = chainBeEnd.getOrigBreakend();

                SvLinkedPair[] matchingPair = {null, null};

                for(int se = SE_START; se <= SE_END; ++se)
                {
                    boolean isStart = isStart(se);
                    SvBreakend chainBe = chain.getOpenBreakend(isStart);

                    if(chainBe == null)
                        continue;

                    chainBe = chainBe.getOrigBreakend();

                    if(chainBe.getSV().getImpliedPloidy() < compDupPloidy * 2)
                        continue;

                    if(getUnlinkedBreakendCount(chainBe) < compDupPloidy)
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

                if(!bestCompDupLinks.isEmpty() && compDupPloidy < maxCompDupPloidy)
                    continue;

                bestCompDupLinks.clear();
                bestCompDupLinks.add(matchingPair[SE_START]);
                bestCompDupLinks.add(matchingPair[SE_END]);
                maxCompDupPloidy = compDupPloidy;
                matchedChain = chain;

                log(LOG_TYPE_VERBOSE, String.format("comDup(%s) ploidy(%d) matched with chain breakends(%s & %s) ploidy(%.1f -> %.1f)",
                        compDup.id(), compDupPloidy, chainBeStart.toString(), chainBeEnd.toString(),
                        chainBeStart.getSV().ploidyMin(), chainBeStart.getSV().ploidyMax()));

                mDiagnostics.logCsv("COMP_DUP", compDup,
                        String.format("ploidy(%.1f-%.1f-%.1f) beStart(%s ploidy=%.1f-%.1f) beEnd(%s ploidy=%.1f-%.1f)",
                                compDup.ploidyMin(), compDup.ploidy(), compDup.ploidyMax(),
                                chainBeStart.toString(), chainBeStart.getSV().ploidyMin(), chainBeStart.getSV().ploidyMax(),
                                chainBeEnd.toString(), chainBeEnd.getSV().ploidyMin(), chainBeEnd.getSV().ploidyMax()));
            }

            // alternatively check if the same SV can be linked to the complex DUP's breakends and satisifies the ploidy constraints
            for(SvLinkedPair pairStart : pairsOnCdStart)
            {
                SvVarData nonFbVar = pairStart.getOtherSV(compDup);
                SvBreakend otherBreakend = pairStart.getOtherBreakend(compDupBeStart);

                if(nonFbVar.getImpliedPloidy() < compDupPloidy * 2)
                    continue;

                // does this exist in the other foldback breakend's set of possible pairs
                for(SvLinkedPair pairEnd : pairsOnCdEnd)
                {
                    SvBreakend otherBreakend2 = pairEnd.getOtherBreakend(compDupBeEnd);

                    // check that available breakends support this SV being connected twice
                    if(otherBreakend.getSV() != otherBreakend2.getSV())
                        continue;

                    if(getUnlinkedBreakendCount(otherBreakend) >= compDupPloidy
                    && getUnlinkedBreakendCount(otherBreakend2) >= compDupPloidy)
                    {
                        if(!bestCompDupLinks.isEmpty() || compDupPloidy > maxCompDupPloidy)
                        {
                            bestCompDupLinks.clear();
                            bestCompDupLinks.add(pairStart);
                            bestCompDupLinks.add(pairEnd);
                            maxCompDupPloidy = compDupPloidy;
                            matchedChain = null;

                            log(LOG_TYPE_VERBOSE, String.format("comDup(%s) ploidy(%d) matched with breakends(%s & %s) ploidy(%.1f -> %.1f)",
                                    compDup.id(), compDupPloidy, otherBreakend, otherBreakend2, nonFbVar.ploidyMin(), nonFbVar.ploidyMax()));
                        }
                    }

                    break;
                }
            }
        }

        if(bestCompDupLinks.isEmpty())
            return Lists.newArrayList();

        return Lists.newArrayList(new ProposedLinks(bestCompDupLinks, maxCompDupPloidy, matchedChain, true));
    }

    private List<ProposedLinks> findFoldbackPairs()
    {
        // both ends of a foldback connect to one end of another SV with ploidy >= 2x

        if(mFoldbacks.isEmpty())
            return Lists.newArrayList();

        List<SvLinkedPair> bestFoldbackLinks = null;
        double maxFoldbackPloidy = 0;
        long foldbackLinkDistance = 0;
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
            double foldbackPloidy = min(getUnlinkedBreakendCount(foldbackStart), getUnlinkedBreakendCount(foldbackEnd));

            if(foldbackPloidy == 0 || foldbackPloidy < maxFoldbackPloidy)
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

                double nonFoldbackPloidy = getUnlinkedBreakendCount(otherBreakend);

                if(nonFoldbackPloidy < foldbackPloidy * 2)
                    continue;

                // does this exist in the other foldback breakend's set of possible pairs
                for(SvLinkedPair pairEnd : pairsOnFbEnd)
                {
                    SvBreakend otherBreakend2 = pairEnd.getOtherBreakend(foldbackEnd);

                    // check that available breakends support this SV being connected twice
                    if(otherBreakend != otherBreakend2)
                        continue;

                    long distance = min(abs(foldbackStart.position() - otherBreakend.position()),
                            abs(foldbackEnd.position() - otherBreakend.position()));

                    if(bestFoldbackLinks.isEmpty() || distance < foldbackLinkDistance)
                    {
                        foldbackLinkDistance = distance;
                        maxFoldbackPloidy = foldbackPloidy;
                        bestFoldbackLinks.clear();
                        bestFoldbackLinks.add(pairStart);
                        bestFoldbackLinks.add(pairEnd);

                        log(LOG_TYPE_VERBOSE, String.format("foldback(%s) ploidy(%d) matched with breakend(%s) ploidy(%.1f -> %.1f)",
                                foldback.id(), foldbackPloidy, otherBreakend, nonFbVar.ploidyMin(), nonFbVar.ploidyMax()));
                    }
                    break;
                }
            }
        }

        if(bestFoldbackLinks.isEmpty())
            return Lists.newArrayList();

        return Lists.newArrayList(new ProposedLinks(bestFoldbackLinks, maxFoldbackPloidy, null, true));
    }

    private List<ProposedLinks> findSingleOptionPairs()
    {
        // find all breakends with only one other link options and order by highest ploidy
        List<ProposedLinks> proposedLinks = Lists.newArrayList();

        for(Map.Entry<SvBreakend, List<SvLinkedPair>> entry : mSvBreakendPossibleLinks.entrySet())
        {
            if(entry.getValue().size() != 1)
                continue;

            SvBreakend limitingBreakend = entry.getKey();

            // special case for DM DUPs - because they can link with themselves at the end, don't restrict their connectivity earlier on
            if(isDoubleMinuteDup() && mDoubleMinuteSVs.contains(limitingBreakend.getSV()))
                continue;

            SvLinkedPair newPair = entry.getValue().get(0);

            if(mSkippedPairs.contains(newPair))
                continue;

            // check for a ploidy overlap

            double minLinkPloidy = getMaxUnlinkedPairCount(newPair);

            if(minLinkPloidy == 0)
                continue;

            // add this to the set of possible if it doesn't clash with any others, and where it does then
            // take the pair with the highest ploidy, following by shortest
            int index = 0;
            boolean canAdd = true;
            while(index < proposedLinks.size())
            {
                final ProposedLinks otherProposedLink = proposedLinks.get(index);
                final SvLinkedPair otherPair = otherProposedLink.Links.get(0);

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

                // one of the breakend is in both the pairs - log a warning and select the highest ploidy option
                SvBreakend otherBreakend = newPair.getOtherBreakend(limitingBreakend);

                if(getUnlinkedBreakendCount(otherBreakend) == 1)
                {
                    log(LOG_TYPE_WARN, String.format("single-option pair(%s len=%d ploidy=%s) clashes with pair(%s len=%d rep=%s)",
                            newPair.toString(), newPair.length(), formatPloidy(minLinkPloidy),
                            otherPair.toString(), otherPair.length(), formatPloidy(otherProposedLink.Ploidy)));
                }

                if(minLinkPloidy > otherProposedLink.Ploidy
                || copyNumbersEqual(minLinkPloidy, otherProposedLink.Ploidy) && newPair.length() < otherPair.length())
                {
                    proposedLinks.remove(index);
                }
                else
                {
                    canAdd = false;
                    break;
                }
            }

            if(canAdd)
            {
                log(LOG_TYPE_VERBOSE, String.format("single-option pair(%s) limited by breakend(%s) minLinkCount(%s)",
                        newPair.toString(), limitingBreakend.toString(), formatPloidy(minLinkPloidy)));

                proposedLinks.add(new ProposedLinks(newPair, minLinkPloidy));
            }
        }

        return proposedLinks;
    }

    private List<ProposedLinks> findAdjacentMatchingPairs()
    {
        if(mAdjacentMatchingPairs.isEmpty())
            return Lists.newArrayList();

        while(!mAdjacentMatchingPairs.isEmpty())
        {
            SvLinkedPair nextPair = mAdjacentMatchingPairs.get(0);

            mAdjacentMatchingPairs.remove(0);

            if(matchesExistingPair(nextPair))
                continue;

            double firstPloidy = getUnlinkedBreakendCount(nextPair.getFirstBreakend());
            double secondPloidy = getUnlinkedBreakendCount(nextPair.getSecondBreakend());

            if(firstPloidy == 0 || secondPloidy == 0)
                continue;

            if(!copyNumbersEqual(firstPloidy, secondPloidy))
                continue;

            // take the average ploidy or calculate a weighted ploidy already?
            // if these links have already been partially used, then incorrect to calculate a weighted ploidy
            double avgPloidy = (firstPloidy + secondPloidy) * 0.5;


            return Lists.newArrayList(new ProposedLinks(nextPair, avgPloidy));
        }

        return Lists.newArrayList();
    }

    private ProposedLinks findMaxReplicationPairs()
    {
        // look at the remaining SVs with replication at least 2 and those with the highest
        // remaining replication count and then those with fewest options
        List<SvLinkedPair> possiblePairs = Lists.newArrayList();

        // first check if there are SVs with a higher replication count, and if so favour these first
        List<SvChainState> maxSvConns = Lists.newArrayList();
        double maxRepCount = 0;

        for(SvChainState svConn : mSvConnectionsMap.values())
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                if (se == SE_END && svConn.SV.isNullBreakend())
                    continue;

                double repCount = svConn.unlinked(se);

                if (repCount <= 1)
                    continue;

                if (repCount > maxRepCount)
                {
                    maxRepCount = repCount;
                    maxSvConns.clear();
                    maxSvConns.add(svConn);
                }
                else if (repCount == maxRepCount)
                {
                    if (!maxSvConns.contains(svConn))
                        maxSvConns.add(svConn);
                }
            }
        }

        if(maxSvConns == null || maxSvConns.isEmpty())
            return null;

        List<SvBreakend> breakendList = Lists.newArrayList();

        for(final SvChainState svConn : maxSvConns)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(se == SE_END && svConn.SV.isNullBreakend())
                    continue;

                if(svConn.unlinked(se) < maxRepCount)
                    continue;

                breakendList.add(svConn.SV.getBreakend(isStart(se)));
            }
        }

        if(breakendList.isEmpty())
            return null;

        if (mLogVerbose)
        {
            for (final SvChainState svConn : maxSvConns)
            {
                LOGGER.debug("restricted to rep SV: {} repCount({})", svConn.SV.id(), svConn.unlinked());
            }
        }

        // next take the pairings with the least alternatives
        // possiblePairs = findFewestOptionPairs(breakendList, true);

        if (possiblePairs.isEmpty())
        {
            // TO-DO prevent this rule from running again

            /*
            // these high-replication SVs yielded no possible links so remove them from consideration
            for (final SvVarData var : maxSvConns)
            {
                log(LOG_TYPE_VERBOSE, String.format("cluster(%s) removing high-replicated SV(%s %s)",
                        mClusterId, var.posId(), var.type()));

                mSvReplicationMap.remove(var);
            }
            */

            // mPairSkipped = true;
        }

        return null;
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

    private ProposedLinks findFewestOptionPairs(List<SvBreakend> breakendList, boolean isRestricted)
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

        // TODO - only select one link, with the highest ploidy or best match
        return new ProposedLinks(minLinkPairs, 0, null, false);
    }

    private void processProposedLinks(List<ProposedLinks> proposedLinksList)
    {
        // now the top candidates to link have been found, take the shortest of them and add this to a chain
        // where possible, add links multiple times according to the min replication of the breakends involved
        // after each link is added, check whether any breakend now has only one link option
        boolean linkAdded = false;

        // boolean isRestrictedSet = mLinkReason == LR_METHOD_ONLY;

        while (!proposedLinksList.isEmpty())
        {
            ProposedLinks proposedLinks = proposedLinksList.get(0);

            proposedLinksList.remove(0);

            linkAdded |= addLinks(proposedLinks);

            if(!mIsValid)
                return;

            /*
            // check whether after adding a link, some SV breakends have only a single possible link
            if(!isRestrictedSet)
            {
                List<SvLinkedPair> restrictedPairs = findSingleOptionPairs();

                // if this finds any conflicting links with the current possible pairs then switch to using the restricted set
                for(SvLinkedPair soPair : restrictedPairs)
                {
                    if(possiblePairs.contains(soPair))
                        continue;

                    for(SvLinkedPair possiblePair : possiblePairs)
                    {
                        if(soPair.hasLinkClash(possiblePair))
                        {
                            possiblePairs = restrictedPairs;
                            isRestrictedSet = true;
                            mLinkReason = LR_METHOD_ONLY;
                            break;
                        }
                    }

                    if(isRestrictedSet)
                        break;
                }
            }

            if(!isRestrictedSet)
            {
                // TODO - doesn't this check for opposite match affect a complex DUP surrounding a DEL?
                // having added a new pair, remove any other conflicting pairs
                int index = 0;
                while (index < possiblePairs.size())
                {
                    SvLinkedPair pair = possiblePairs.get(index);
                    if (!mHasReplication && pair.oppositeMatch(shortestPair))
                    {
                        possiblePairs.remove(index);
                        continue;
                    }

                    if (mHasReplication)
                    {
                        if(getUnlinkedBreakendCount(pair.getBreakend(true)) == 0
                        || getUnlinkedBreakendCount(pair.getBreakend(false)) == 0)
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
            */
        }

        if(linkAdded)
        {
            mSkippedPairs.clear(); // any skipped links can now be re-evaluated
        }
    }

    private static int SPEC_LINK_INDEX = -1;
    // private static int SPEC_LINK_INDEX = 26;

    private boolean addLinks(final ProposedLinks proposedLinks)
    {
        // if a chain is specified, add the links to it
        // otherwise look for a chain which can link in these new pairs
        // and if none can be found, create a new chain with them

        // if no chain has a ploidy matching that of the new link and the new link is lower, then split the chain
        // a chain has a matching ploidy then recalculate it with the new SV's ploidy and uncertainty

        SvLinkedPair newPair = proposedLinks.Links.get(0);

        boolean[] pairLinkedOnFirst = {false, false};
        boolean addedToChain = false;
        boolean addToStart = false;
        boolean linkClosesChain = false;

        SvChain targetChain = null;

        if(proposedLinks.ChainTarget != null)
        {
            targetChain = proposedLinks.ChainTarget;
        }
        else
        {
            for (SvChain chain : mChains)
            {
                boolean canAddToStart = chain.canAddLinkedPair(newPair, true, true);
                boolean canAddToEnd = chain.canAddLinkedPair(newPair, false, true);

                if (!canAddToStart && !canAddToEnd)
                    continue;

                boolean couldCloseChain = canAddToStart && canAddToEnd;

                if (couldCloseChain)
                {
                    if (isDoubleMinuteDup() && mSvConnectionsMap.size() == 1 && mSvConnectionsMap.get(mDoubleMinuteSVs.get(0)) != null)
                    {
                        // allow the chain to be closed if this is the last pair other than excess DM DUP replicated SVs

                    }
                    else
                    {
                        // the link can be added to both ends, which would close the chain - so search for an alternative SV on either end
                        // to keep it open while still adding the link

                        log(LOG_TYPE_VERBOSE, String.format("skipping linked pair(%s) would close existing chain(%d)",
                                newPair.toString(), chain.id()));

                        if (!mSkippedPairs.contains(newPair))
                        {
                            mPairSkipped = true;
                            mSkippedPairs.add(newPair);
                        }

                        linkClosesChain = true;
                        continue;
                    }
                }

                boolean linkOnFirst = newPair.first() == chain.getFirstSV() || newPair.first() == chain.getLastSV();

                // test ploidy match
                boolean ploidyMatched = false;

                if ((linkOnFirst && copyNumbersEqual(chain.ploidy(), newPair.first().ploidy()))
                    || (!linkOnFirst && copyNumbersEqual(chain.ploidy(), newPair.second().ploidy())))
                {
                    ploidyMatched = true;
                }
                else
                {
                    // this may be a candidate for splitting a chain

                }
            }
        }

        if(targetChain != null)
        {

            targetChain.addLink(newPair, addToStart);
            addedToChain = true;

            LOGGER.debug("index({}) method({}) adding linked pair({} {} len={}) to existing chain({}) {}",
                    mLinkIndex, mLinkReason, newPair.toString(), newPair.assemblyInferredStr(), newPair.length(),
                    targetChain.id(), addToStart ? "start" : "end");
        }
        else
        {
            if(linkClosesChain)
                return false; // skip this link for now

            SvChain chain = new SvChain(mNextChainId++);
            mChains.add(chain);

            chain.addLink(newPair, true);

            PloidyCalcData ploidyData = calcPloidyUncertainty(new PloidyCalcData(newPair.first()), new PloidyCalcData(newPair.second()));
            chain.setPloidyData(ploidyData.PloidyEstimate, ploidyData.PloidyUncertainty);

            pairLinkedOnFirst[SE_START] = true;
            pairLinkedOnFirst[SE_END] = true;

            LOGGER.debug("index({}) method({}) adding linked pair({} {}) to new chain({}) new ploidy({})",
                    mLinkIndex, mLinkReason, newPair.toString(), newPair.assemblyInferredStr(), chain.id(),
                    String.format("%.1f unc=%.1f", ploidyData.PloidyEstimate, ploidyData.PloidyUncertainty));
        }

        newPair.setLinkReason(mLinkReason, mLinkIndex);

        registerNewLink(newPair, pairLinkedOnFirst, proposedLinks.Ploidy);
        ++mLinkIndex;

        if(mRunValidation)
            mDiagnostics.checkHasValidState(mLinkIndex);

        if(addedToChain)
        {
            // now see if any partial chains can be linked
            reconcileChains();
        }

        return true;
    }

    private void registerNewLink(final SvLinkedPair newPair, boolean[] pairToChain, double linkPloidy)
    {
        for (int be = SE_START; be <= SE_END; ++be)
        {
            boolean isStart = isStart(be);

            final SvBreakend breakend = newPair.getBreakend(isStart);
            final SvBreakend otherPairBreakend = newPair.getOtherBreakend(breakend);
            final SvVarData var = breakend.getSV();

            SvChainState svConn = mSvConnectionsMap.get(var);

            if(svConn == null || svConn.maxUnlinked(breakend.usesStart()) < MIN_PLOIDY)
            {
                LOGGER.error("breakend({}) connections exhausted: {}",
                        breakend.toString(), svConn != null ? svConn.toString() : "null");
                mIsValid = false;
                return;
            }

            svConn.add(breakend.usesStart(), linkPloidy);
            svConn.addConnection(otherPairBreakend, breakend.usesStart());

            final SvBreakend otherSvBreakend = var.getBreakend(!breakend.usesStart());

            boolean hasUnlinkedBreakend = true;

            boolean svExhausted = false;
            if(svConn.maxUnlinked(breakend.usesStart()) == 0)
            {
                // mUnlinkedBreakendMap.remove(origBreakend);
                hasUnlinkedBreakend = false;

                // LOGGER.debug("breakend({}) has no more possible links", breakend);

                if(svConn.maxUnlinked(!breakend.usesStart()) == 0)
                {
                    svExhausted = true;

                    if(var.isFoldback())
                    {
                        // remove if no other instances of this SV remain
                        mFoldbacks.remove(var);
                    }
                    else if(mComplexDupCandidates.contains(var))
                    {
                        mComplexDupCandidates.remove(var);
                    }
                }
            }

            List<SvLinkedPair> possibleLinks = mSvBreakendPossibleLinks.get(breakend);

            if (possibleLinks != null && !hasUnlinkedBreakend)
            {
                // not more replicated breakends exist so any possible links that depend on it
                removePossibleLinks(possibleLinks, breakend);
            }

            if(!mHasReplication)
            {
                // check for an opposite pairing between these 2 SVs - need to look into other breakends' lists
                possibleLinks = otherSvBreakend != null && !mSvBreakendPossibleLinks.isEmpty()
                        ? mSvBreakendPossibleLinks.get(otherSvBreakend)
                        : null;

                if (possibleLinks != null)
                {
                    final SvBreakend otherOrigBreakendAlt = otherPairBreakend.getOtherBreakend();

                    if (otherOrigBreakendAlt != null)
                    {
                        for (SvLinkedPair pair : possibleLinks)
                        {
                            if (pair.hasBreakend(otherSvBreakend) && pair.hasBreakend(otherOrigBreakendAlt))
                            {
                                possibleLinks.remove(pair);
                                break;
                            }
                        }
                    }
                }
            }

            if(mHasReplication)
            {
                // reduce replication counts for breakends which are added to a chain
                if (pairToChain[be])
                {
                    checkSvComplete(svConn);
                }
            }
        }

        // track unique pairs to avoid conflicts (eg end-to-end and start-to-start)
        if(!matchesExistingPair(newPair))
        {
            mUniquePairs.add(newPair);
        }
    }

    private static final double MIN_PLOIDY = 0.2;

    private void checkSvComplete(final SvChainState svConn)
    {
        if(svConn.maxUnlinked(true) < MIN_PLOIDY && (svConn.SV.isNullBreakend() || svConn.maxUnlinked(false) < MIN_PLOIDY))
        {
            log(LOG_TYPE_VERBOSE, String.format("SV(%s) connections exhausted", svConn.toString()));
            mSvConnectionsMap.remove(svConn.SV);
            mSvCompletedConnections.add(svConn);
        }
    }

    /*
    private SvBreakend findUnlinkedMatchingBreakend(final SvBreakend breakend)
    {
        // get the next available breakend (thereby reducing the replicated instances)
        final List<SvBreakend> breakendList = mUnlinkedBreakendMap.get(breakend.getOrigBreakend());

        if(breakendList == null || breakendList.isEmpty())
            return null;

        return breakendList.get(0);
    }
    */

    private double getUnlinkedBreakendCount(final SvBreakend breakend)
    {
        // List<SvBreakend> beList = mUnlinkedBreakendMap.get(breakend);
        // return beList != null ? beList.size() : 0;
        SvChainState svConn = mSvConnectionsMap.get(breakend.getSV());
        if(svConn == null)
            return 0;

        double breakendPloidy = svConn.unlinked(breakend.usesStart());

        return breakendPloidy >= MIN_PLOIDY ? breakendPloidy : 0;
    }

    private double getMaxUnlinkedPairCount(final SvLinkedPair pair)
    {
        double first = getUnlinkedBreakendCount(pair.getBreakend(true));
        double second = getUnlinkedBreakendCount(pair.getBreakend(false));
        return min(first, second);
    }

    private double getMaxUnlinkedPairsCount(final SvLinkedPair pair1, final SvLinkedPair pair2)
    {
        // checks for repeated breakends
        List<SvBreakend> uniqueBreakends = Lists.newArrayList();
        double minCount = 0;

        for(int i = 0; i <= 1; ++i)
        {
            SvLinkedPair pair = (i == 0) ? pair1 : pair2;

            for (int be = SE_START; be <= SE_END; ++be)
            {
                SvBreakend breakend = pair.getBreakend(isStart(be));
                double count = getUnlinkedBreakendCount(breakend);

                if (count == 0)
                    return 0;

                if (uniqueBreakends.isEmpty())
                {
                    minCount = count;
                    uniqueBreakends.add(breakend);
                    continue;
                }

                if (uniqueBreakends.contains(breakend))
                {
                    count = count / 2;
                }
                else
                {
                    uniqueBreakends.add(breakend);
                }

                minCount = min(minCount, count);
            }
        }

        return minCount;
    }

    private void addAssemblyLinksToChains()
    {
        if(mAssembledLinks.isEmpty())
            return;

        mLinkReason = LR_METHOD_ASMB;

        for(SvLinkedPair pair : mAssembledLinks)
        {
            if(!mHasReplication)
            {
                double minPloidy = getMaxUnlinkedPairCount(pair);
                addLinks(new ProposedLinks(pair, minPloidy));
                continue;
            }

            // replicate any assembly links where the ploidy supports it, taking note of multiple connections between the same
            // breakend and other breakends eg if a SV has ploidy 2 and 2 different assembly links, it can only link once, whereas
            // if it has ploidy 2 and 1 link it should be made twice, and any higher combinations are unclear
            double[] repCounts = new double[SE_PAIR];
            boolean[] hasOtherMultiPloidyLinks = {false, false};
            boolean[] hasOtherMultiAssemblyLinks = {false, false};
            int[] assemblyLinkCount = new int[SE_PAIR];
            double pairPloidy = 0;

            for (int be = SE_START; be <= SE_END; ++be)
            {
                boolean isStart = isStart(be);

                final SvBreakend breakend = pair.getBreakend(isStart);
                repCounts[be] = getUnlinkedBreakendCount(breakend);

                final List<SvLinkedPair> assemblyLinks = breakend.getSV().getAssembledLinkedPairs(breakend.usesStart());
                assemblyLinkCount[be] = assemblyLinks.size();

                if(assemblyLinkCount[be] > 1)
                {
                    for(final SvLinkedPair assemblyLink : assemblyLinks)
                    {
                        final SvBreakend otherBreakend = assemblyLink.getOtherBreakend(breakend);

                        if(getUnlinkedBreakendCount(otherBreakend) > 1)
                            hasOtherMultiPloidyLinks[be] = true;

                        if(otherBreakend.getSV().getAssembledLinkedPairs(otherBreakend.usesStart()).size() > 1)
                            hasOtherMultiAssemblyLinks[be] = true;
                    }
                }
            }

            // most likely scenario first
            if(assemblyLinkCount[SE_START] == 1 && assemblyLinkCount[SE_END] == 1)
            {
                pairPloidy = min(repCounts[SE_START], repCounts[SE_END]);
            }
            else if(repCounts[SE_START] >= 2 && repCounts[SE_END] >= 2
                && (!hasOtherMultiPloidyLinks[SE_START] || !hasOtherMultiAssemblyLinks[SE_START])
                && (!hasOtherMultiPloidyLinks[SE_END] || !hasOtherMultiAssemblyLinks[SE_END]))
            {
                // both SVs allow for 2 or more repeats, so if the max other assembled link SV ploidies are all <= 1,
                // then these links can safely be repeated
                int firstOtherLinks = assemblyLinkCount[SE_START] - 1;
                int secondOtherLinks = assemblyLinkCount[SE_END] - 1;

                pairPloidy = min(repCounts[SE_START] - firstOtherLinks, repCounts[SE_END] - secondOtherLinks);
            }

            if(pairPloidy > 1)
            {
                LOGGER.debug("assembly pair({}) ploidy({}): first(rep={} links={}) second(rep={} links={})",
                        pair.toString(), formatPloidy(pairPloidy), repCounts[SE_START], assemblyLinkCount[SE_START],
                        repCounts[SE_END], assemblyLinkCount[SE_END]);
            }

            addLinks(new ProposedLinks(pair, pairPloidy));
        }

        if(!mChains.isEmpty())
        {
            LOGGER.debug("created {} partial chains from {} assembly links", mChains.size(), mAssembledLinks.size());
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
                                mSvBreakendPossibleLinks.remove(otherBreakend);
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
            mSvBreakendPossibleLinks.remove(origBreakend);
        }
    }

    private int getClusterChrBreakendIndex(final SvBreakend breakend)
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

        for (final Map.Entry<String, List<SvBreakend>> entry : mChrBreakendMap.entrySet())
        {
            final String chromosome = entry.getKey();
            final List<SvBreakend> breakendList = entry.getValue();
            final double[][] allelePloidies = mClusterPloidyLimits.getChrAllelePloidies().get(chromosome);

            for (int i = 0; i < breakendList.size() -1; ++i)
            {
                final SvBreakend lowerBreakend = breakendList.get(i);

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

                double lowerPloidy = getUnlinkedBreakendCount(lowerBreakend);

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

                    // make note of any pairs formed from adjacent facing breakends
                    if(j == i + 1 && copyNumbersEqual(lowerPloidy, getUnlinkedBreakendCount(upperBreakend)))
                    {
                        mAdjacentMatchingPairs.add(newPair);
                    }

                    lowerPairs.add(newPair);

                    List<SvLinkedPair> upperPairs = mSvBreakendPossibleLinks.get(upperBreakend);

                    if(upperPairs == null)
                    {
                        upperPairs = Lists.newArrayList();
                        mSvBreakendPossibleLinks.put(upperBreakend, upperPairs);
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
                            log(LOG_TYPE_VERBOSE, String.format("breakends lower(%d: %s) limited at upper(%d: %s) with clusterAP(%.2f)",
                                    i, lowerBreakend.toString(), j, upperBreakend.toString(), clusterAP));

                            break;
                        }
                    }
                }
            }
        }
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
                    if(otherSV == higherPloidyBreakend.getSV())
                    {
                        log(LOG_TYPE_INFO, String.format("identified complex dup(%s %s) ploidy(%.1f -> %.1f) vs SV(%s) ploidy(%.1f -> %.1f)",
                                var.posId(), var.type(), var.ploidyMin(), var.ploidyMax(), higherPloidyBreakend.getSV().id(),
                                higherPloidyBreakend.getSV().ploidyMin(), higherPloidyBreakend.getSV().ploidyMax()));
                    }
                    else
                    {
                        log(LOG_TYPE_INFO, String.format("identified complex dup(%s %s) ploidy(%.1f -> %.1f) vs SV(%s) ploidy(%.1f -> %.1f) & SV(%s) ploidy(%.1f -> %.1f)",
                                var.posId(), var.type(), var.ploidyMin(), var.ploidyMax(),
                                otherSV.id(), otherSV.ploidyMin(), otherSV.ploidyMax(), higherPloidyBreakend.getSV().id(),
                                higherPloidyBreakend.getSV().ploidyMin(), higherPloidyBreakend.getSV().ploidyMax()));
                    }

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

    private void populateSvPloidyMap()
    {
        // make a cache of all unchained breakends in those of replicated SVs
        for(final SvVarData var : mSvList)
        {
            mSvConnectionsMap.put(var, new SvChainState(var, !mHasReplication));
        }
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

    private boolean isDoubleMinuteDup()
    {
        return mDoubleMinuteSVs.size() == 1 && mDoubleMinuteSVs.get(0).type() == DUP;
    }

    private void checkDoubleMinuteChains()
    {
        // if there is a single chain which contains all DM SVs, attempt to close the chain
        if(mDoubleMinuteSVs.isEmpty())
            return;

        if(mChains.size() != 1)
            return;

        SvChain chain = mChains.get(0);

        // allow any excess breakends from a single DM DUP to be added to the chain
        /* probably no point since just creates replicated identical links and won't impact visualisation

        if(isDoubleMinuteDup())
        {
            final SvVarData dupDM = mDoubleMinuteSVs.get(0);
            final SvBreakend dupStart = dupDM.getBreakend(true);
            final SvBreakend dupEnd = dupDM.getBreakend(false);

            int remainingBreakends = min(getUnlinkedBreakendCount(dupStart), getUnlinkedBreakendCount(dupEnd));

            if(remainingBreakends > 0)
            {
                LOGGER.debug("cluster({}) adding DUP pair to DM chain {} times", mClusterId, remainingBreakends);

                if(chain.getFirstSV().equals(dupDM, true) || chain.getLastSV().equals(dupDM, true))
                {
                    final List<SvBreakend> startBreakendList = mUnlinkedBreakendMap.get(dupStart);
                    final List<SvBreakend> endBreakendList = mUnlinkedBreakendMap.get(dupEnd);

                    // work out which end of the chain has this DUP if any
                    boolean linkOnStart = chain.getFirstSV().equals(dupDM, true);

                    for(int i = 0; i < remainingBreakends; ++i)
                    {
                        SvBreakend chainBreakend = chain.getOpenBreakend(linkOnStart);

                        SvBreakend otherBreakendOrig = dupStart == chainBreakend.getOrigBreakend() ? dupEnd : dupStart;
                        SvBreakend otherBreakend = findUnlinkedMatchingBreakend(otherBreakendOrig);

                        if(otherBreakend == null || chainBreakend.getSV() == otherBreakend.getSV())
                            break;

                        SvLinkedPair newLink = SvLinkedPair.from(chainBreakend, otherBreakend);

                        chain.addLink(newLink, linkOnStart);
                        newLink.setLinkReason(LR_METHOD_DM_DUP, mLinkIndex++);

                        if(chainBreakend.usesStart())
                        {
                            startBreakendList.remove(chainBreakend);
                            endBreakendList.remove(otherBreakend);
                        }
                        else
                        {
                            startBreakendList.remove(otherBreakend);
                            endBreakendList.remove(chainBreakend);
                        }
                    }
                }
                else
                {
                    // just add the extra links even though they're not in the correct location
                    for(int i = 0; i < remainingBreakends; ++i)
                    {
                        SvLinkedPair newLink = SvLinkedPair.from(dupStart, dupEnd);

                        chain.addLink(newLink, 0);
                        newLink.setLinkReason(LR_METHOD_DM_DUP, mLinkIndex++);
                    }
                }
            }
        }
        */

        int chainedDmSVs = (int)mDoubleMinuteSVs.stream().filter(x -> chain.hasSV(x, true)).count();

        if(chainedDmSVs != mDoubleMinuteSVs.size())
            return;

        SvBreakend chainStart = chain.getOpenBreakend(true);
        SvBreakend chainEnd = chain.getOpenBreakend(false);

        if(chainStart != null && !chainStart.getSV().isNullBreakend() && chainEnd != null && !chainEnd.getSV().isNullBreakend())
        {
            if (areLinkedSection(chainStart.getSV(), chainEnd.getSV(), chainStart.usesStart(), chainEnd.usesStart(), false))
            {
                SvLinkedPair pair = SvLinkedPair.from(chainStart, chainEnd);

                if (chain.linkWouldCloseChain(pair))
                {
                    chain.addLink(pair, true);
                    pair.setLinkReason(LR_METHOD_DM_CLOSE, mLinkIndex++);

                    LOGGER.debug("cluster({}) closed DM chain", mClusterId);
                }
            }
        }
    }

    private void log(int level, final String message)
    {
        if(level >= LOG_TYPE_VERBOSE && !mLogVerbose)
            return;

        if(level >= LOG_TYPE_INFO && !LOGGER.isDebugEnabled())
            return;

        if(level <= LOG_TYPE_INFO)
        {
            mDiagnostics.addMessage(level, message);
        }

        LOGGER.debug(message);
    }

}
