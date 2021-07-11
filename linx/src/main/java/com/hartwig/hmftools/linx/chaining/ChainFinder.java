package com.hartwig.hmftools.linx.chaining;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INF;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.LINE_CHAINS;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.runAnnotation;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;
import static com.hartwig.hmftools.linx.chaining.ChainJcnLimits.CLUSTER_ALLELE_JCN_MIN;
import static com.hartwig.hmftools.linx.chaining.ChainJcnLimits.jcnExceedsMajorAlleleJcn;
import static com.hartwig.hmftools.linx.chaining.ChainJcnLimits.jcnMatchForSplits;
import static com.hartwig.hmftools.linx.chaining.ChainLinkAllocator.belowJcnThreshold;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.identicalChain;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.reconcileChains;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.getMinTemplatedInsertionLength;
import static com.hartwig.hmftools.linx.chaining.LinkSkipType.JCN_MISMATCH;
import static com.hartwig.hmftools.linx.chaining.SvChain.checkIsValid;
import static com.hartwig.hmftools.linx.types.ResolvedType.LINE;

import static org.apache.logging.log4j.Level.TRACE;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;

/* ChainFinder - forms one or more chains from the SVs in a cluster

    Set-up:
    - form all assembled links and connect into chains
    - identify foldbacks and complex DUP-type SVs, since these will be prioritised during chaining
    - create a cache of all possible linked pairs

    Routine:
    - apply priority rules to find the next possible link(s)
    - add the link to an existing chain or a new chain if required
    - remove the breakends & link from further consideration
    - repeat until no further links can be made

    Priority rules in order:
    - Foldbacks - look for a breakend which can link to both ends of a foldback
    - Single-Option - if a breakend has only one possible link, select this one
    - Ploidy-Match - starting with the highest ploidy SV, only link SVs of the same ploidy
    - Shortest - after all other rules, if there is more than 1 possible link then choose the shortest
*/

public class ChainFinder
{
    private int mClusterId;
    public String mSampleId;
    private boolean mHasReplication;
    private boolean mIsLineCluster;

    // input state
    private final List<SvVarData> mSvList;
    private final List<SvVarData> mFoldbacks;
    private final List<SvVarData> mDoubleMinuteSVs;
    private final List<LinkedPair> mAssembledLinks;
    private Map<String,List<SvBreakend>> mChrBreakendMap;

    // links a breakend to its position in a chromosome's breakend list - only applicable if a subset of SVs are being chained
    private boolean mIsClusterSubset;
    private final Map<SvBreakend,Integer> mSubsetBreakendClusterIndexMap;

    // chaining state
    private final Map<SvVarData,List<LinkedPair>> mComplexDupCandidates; // identified SVs which duplication another SV
    private final List<LinkedPair> mAdjacentMatchingPairs;
    private final List<LinkedPair> mAdjacentPairs;

    private final List<SvChain> mChains;
    private final List<SvChain> mUniqueChains;

    private final ChainRuleSelector mRuleSelector;
    private final ChainLinkAllocator mLinkAllocator;
    private final LineChainer mLineChainer;

    // a cache of cluster ploidy boundaries which links cannot cross
    private final ChainJcnLimits mClusterJcnLimits;

    // determined up-front - the set of all possible links from a specific breakend to other breakends, ordered shortest first
    private final Map<SvBreakend, List<LinkedPair>> mSvBreakendPossibleLinks;

    private List<SvVarData> mReplicatedSVs;
    private List<SvBreakend> mReplicatedBreakends;

    public static final double MIN_CHAINING_JCN_LEVEL = 0.05;

    private boolean mIsValid;
    private boolean mLogVerbose;
    private Level mLogLevel;
    private boolean mRunValidation;
    private boolean mUseAlleleJCNs;

    public static final String LR_METHOD_DM_CLOSE = "DM_CLOSE";

    // self-analysis only
    private final ChainDiagnostics mDiagnostics;

    public ChainFinder()
    {
        mSvList = Lists.newArrayList();
        mFoldbacks = Lists.newArrayList();
        mDoubleMinuteSVs = Lists.newArrayList();
        mIsClusterSubset = false;
        mAssembledLinks = Lists.newArrayList();
        mChrBreakendMap = null;

        mClusterJcnLimits = new ChainJcnLimits();

        mAdjacentMatchingPairs = Lists.newArrayList();
        mAdjacentPairs = Lists.newArrayList();
        mSubsetBreakendClusterIndexMap = Maps.newHashMap();
        mComplexDupCandidates = Maps.newHashMap();
        mChains = Lists.newArrayList();
        mUniqueChains = Lists.newArrayList();
        mSvBreakendPossibleLinks = Maps.newHashMap();
        mReplicatedSVs = Lists.newArrayList();
        mReplicatedBreakends = Lists.newArrayList();

        mLinkAllocator = new ChainLinkAllocator(mClusterJcnLimits, mSvBreakendPossibleLinks, mChains, mDoubleMinuteSVs);

        mRuleSelector = new ChainRuleSelector(mLinkAllocator, mClusterJcnLimits,
                mSvBreakendPossibleLinks, mFoldbacks, mComplexDupCandidates,
                mAdjacentMatchingPairs, mAdjacentPairs, mChains);

        mLineChainer = new LineChainer();

        mHasReplication = false;
        mLogVerbose = false;
        mLogLevel = LNX_LOGGER.getLevel();
        mRunValidation = false;
        mIsValid = true;
        mSampleId= "";
        mUseAlleleJCNs = false;

        mDiagnostics = new ChainDiagnostics(
                mLinkAllocator.getSvConnections(), mLinkAllocator.getSvCompletedConnections(), mChains, mUniqueChains,
                mSvBreakendPossibleLinks, mDoubleMinuteSVs, mLinkAllocator.getUniquePairs());
    }

    public void initialiseOutput(final LinxConfig config)
    {
        mDiagnostics.setOutputDir(config.OutputDataPath, config.Output.LogChainingMaxSize);

        if(runAnnotation(config.RequiredAnnotations, LINE_CHAINS))
            mLineChainer.initialiseOutput(config.OutputDataPath);
    }

    public void clear()
    {
        mClusterId = -1;
        mIsLineCluster = false;
        mSvList.clear();
        mFoldbacks.clear();
        mDoubleMinuteSVs.clear();
        mIsClusterSubset = false;
        mAssembledLinks.clear();
        mChrBreakendMap = null;

        mAdjacentMatchingPairs.clear();
        mAdjacentPairs.clear();
        mSubsetBreakendClusterIndexMap.clear();
        mComplexDupCandidates.clear();
        mChains.clear();
        mUniqueChains.clear();
        mSvBreakendPossibleLinks.clear();
        mReplicatedSVs.clear();
        mReplicatedBreakends.clear();

        mIsValid = true;

        mDiagnostics.clear();
        mLineChainer.clear();
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

        mClusterId = cluster.id();
        mIsLineCluster = cluster.getResolvedType() == LINE;

        if(mIsLineCluster)
            mLineChainer.initialise(mSampleId, cluster);

        mSvList.addAll(cluster.getSVs());
        mFoldbacks.addAll(cluster.getFoldbacks());
        mDoubleMinuteSVs.addAll(cluster.getDoubleMinuteSVs());
        mAssembledLinks.addAll(cluster.getAssemblyLinkedPairs());
        mChrBreakendMap = cluster.getChrBreakendMap();
        mHasReplication = cluster.requiresReplication();
        mIsClusterSubset = false;
        mLinkAllocator.initialise(mClusterId);
        mRuleSelector.initialise(mClusterId, mHasReplication);
    }

    public void initialise(SvCluster cluster, final List<SvVarData> svList, boolean requiresReplication)
    {
        // chain a specific subset of a cluster's SVs - current used to test double-minute completeness
        clear();

        mIsClusterSubset = true;
        mClusterId = cluster.id();
        mIsLineCluster = false;

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

        mHasReplication = requiresReplication;

        for(SvVarData var : svList)
        {
            mSvList.add(var);

            if(var.isFoldback() && mFoldbacks.contains(var))
                mFoldbacks.add(var);

            for(int se = SE_START; se <= SE_END; ++se)
            {
                // only add an assembled link if it has a partner in the provided SV set, and can be replicated equally
                for (LinkedPair link : var.getAssembledLinkedPairs(isStart(se)))
                {
                    final SvVarData otherVar = link.getOtherSV(var);

                    if(!svList.contains(otherVar))
                        continue;

                    mAssembledLinks.add(link);
                }
            }
        }

        mLinkAllocator.initialise(mClusterId);
        mRuleSelector.initialise(mClusterId, mHasReplication);

        mDoubleMinuteSVs.addAll(svList);
    }

    public void setRunValidation(boolean toggle) { mRunValidation = toggle; }
    public void setUseAllelePloidies(boolean toggle) { mUseAlleleJCNs = toggle; }

    public final List<SvChain> getUniqueChains()
    {
        return mUniqueChains;
    }
    public double getValidAllelePloidySegmentPerc() { return mClusterJcnLimits.getValidAlleleJcnSegmentPerc(); }
    public final long[] calcRangeData() { return mClusterJcnLimits.calcRangeData(); }
    public final ChainDiagnostics getDiagnostics() { return mDiagnostics; }

    public void close()
    {
        mDiagnostics.close();
        mLineChainer.close();
    }

    public void formChains(boolean assembledLinksOnly)
    {
        if(mIsLineCluster)
        {
            formLineChains();
            return;
        }

        if(mSvList.size() < 2)
            return;

        if(mSvList.size() >= 4)
        {
            LNX_LOGGER.debug("cluster({}) starting chaining with assemblyLinks({}) svCount({}) replication({})",
                    mClusterId, mAssembledLinks.size(), mSvList.size(), mHasReplication);
        }

        enableLogVerbose();

        mClusterJcnLimits.initialise(mClusterId, mChrBreakendMap);

        buildChains(assembledLinksOnly);

        checkChains();
        removeIdenticalChains();

        mDiagnostics.chainingComplete();

        disableLogVerbose();

        if(!isValid())
        {
            LNX_LOGGER.warn("cluster({}) chain finding failed", mClusterId);
        }
    }

    private void formLineChains()
    {
        mLineChainer.formChains();
        mUniqueChains.addAll(mLineChainer.getChains());
    }

    private boolean isValid()
    {
        return mIsValid && mLinkAllocator.isValid();
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
                if (identicalChain(chain, newChain, false))
                {
                    LNX_LOGGER.debug("cluster({}) skipping duplicate chain({}) ploidy({}) vs origChain({}) ploidy({})",
                            mClusterId, newChain.id(), formatJcn(newChain.jcn()), chain.id(), formatJcn(chain.jcn()));

                    // combine the ploidies
                    chain.setJcnData(chain.jcn() + newChain.jcn(), chain.jcnUncertainty());
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
        // add these chains to the cluster, but skip any which are identical to existing ones,
        // which can happen for clusters with replicated SVs
        mUniqueChains.stream().forEach(chain -> checkAddNewChain(chain, cluster));

        for(int i = 0; i < cluster.getChains().size(); ++i)
        {
            final SvChain chain = cluster.getChains().get(i);

            if(LNX_LOGGER.isDebugEnabled())
            {
                LNX_LOGGER.debug("cluster({}) added chain({}) ploidy({}) with {} linked pairs:",
                        mClusterId, chain.id(), formatJcn(chain.jcn()), chain.getLinkCount());
                chain.logLinks();
            }

            chain.setId(i); // set after logging so can compare with logging during building
        }
    }

    private void checkAddNewChain(final SvChain newChain, SvCluster cluster)
    {
        if(!mHasReplication)
        {
            cluster.addChain(newChain, false);
            return;
        }

        // any identical chains (including precise subsets of longer chains) will have their replicated SVs entirely removed
        if(!mUniqueChains.contains(newChain))
            return;

        cluster.addChain(newChain, false);
    }

    private void buildChains(boolean assembledLinksOnly)
    {
        mLinkAllocator.populateSvJcnMap(mSvList, mHasReplication);

        mDiagnostics.initialise(mClusterId, mHasReplication);

        // first make chains out of any assembly links
        mLinkAllocator.addAssemblyLinksToChains(mAssembledLinks, mHasReplication);

        if(assembledLinksOnly)
            return;

        if(mUseAlleleJCNs && mHasReplication)
            mClusterJcnLimits.determineBreakendJCNs();

        determinePossibleLinks();

        mDiagnostics.setPriorityData(Lists.newArrayList(mComplexDupCandidates.keySet()), mFoldbacks);

        int iterationsWithoutNewLinks = 0; // protection against loops

        while (true)
        {
            mLinkAllocator.clearSkippedState();
            int lastAddedIndex = mLinkAllocator.getLinkIndex();

            List<ProposedLinks> proposedLinks = mRuleSelector.findProposedLinks();

            if(proposedLinks.isEmpty())
            {
                if(!mLinkAllocator.pairSkipped())
                    break;
            }
            else
            {
                mLinkAllocator.processProposedLinks(proposedLinks);

                if (mRunValidation)
                {
                    checkChains();
                    mDiagnostics.checkHasValidState(mLinkAllocator.getLinkIndex());
                }

                if(!isValid())
                    return;
            }

            // as a safety check, exit if no link is allocated (due to skipping) from too many attempts
            if(lastAddedIndex == mLinkAllocator.getLinkIndex())
            {
                ++iterationsWithoutNewLinks;

                if (iterationsWithoutNewLinks >= 50)
                {
                    LNX_LOGGER.info("cluster({}) {} iterations without adding a link, skipped pairs: closing({}) ploidyMismatch({})",
                            mClusterId, iterationsWithoutNewLinks,
                            mLinkAllocator.getSkippedPairCount(LinkSkipType.CLOSING),
                            mLinkAllocator.getSkippedPairCount(LinkSkipType.JCN_MISMATCH));

                    // mIsValid = false;
                    break;
                }
            }
            else
            {
                iterationsWithoutNewLinks = 0;
            }

            mDiagnostics.checkProgress(mLinkAllocator.getLinkIndex());
        }

        // TEMP - analysis of links not made due to chain ploidy mismatches where the SVs themselves would satisfy the same test
        if(mLinkAllocator.getSkippedPairCount(LinkSkipType.JCN_MISMATCH) > 0)
        {
            int ploidyMismatches = 0;
            int ploidyOverlaps = 0;
            for(Map.Entry<LinkedPair,LinkSkipType> entry : mLinkAllocator.getSkippedPairs().entrySet())
            {
                if(entry.getValue() != JCN_MISMATCH)
                    continue;

                ++ploidyMismatches;

                final LinkedPair pair = entry.getKey();

                if(ChainJcnLimits.jcnMatch(pair.firstBreakend(), pair.secondBreakend()))
                {
                    ++ploidyOverlaps;
                    LNX_LOGGER.debug("skipped pair({}) has ploidy overlap be1({} & {}) and be2({} & {})",
                            pair.toString(), formatJcn(pair.firstBreakend().jcn()), formatJcn(pair.firstBreakend().jcnUncertainty()),
                            formatJcn(pair.secondBreakend().jcn()), formatJcn(pair.secondBreakend().jcnUncertainty()));
                }
            }

            if(ploidyOverlaps > 0)
            {
                LNX_LOGGER.info("sample({}) cluster({}) ploidy skips({}) with overlaps({}) chainCount({})",
                        mSampleId, mClusterId, ploidyMismatches, ploidyOverlaps, mChains.size());
            }
        }

        if(mChains.size() < 50)
        {
            reconcileChains(mChains, false, mLinkAllocator.getNextChainId(), true);
        }

        checkDoubleMinuteChains();
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
        // do not chain past a zero cluster allele JCN
        // identify potential complex DUP candidates along the way

        for (final Map.Entry<String, List<SvBreakend>> entry : mChrBreakendMap.entrySet())
        {
            final String chromosome = entry.getKey();
            final List<SvBreakend> breakendList = entry.getValue();
            final List<SegmentJcn> alleleJCNs = mClusterJcnLimits.getChrAlleleJCNs().get(chromosome);

            for (int i = 0; i < breakendList.size() -1; ++i)
            {
                final SvBreakend lowerBreakend = breakendList.get(i);

                if(belowJcnThreshold(lowerBreakend.getSV()))
                    continue;

                if(lowerBreakend.orientation() != -1)
                    continue;

                if(alreadyLinkedBreakend(lowerBreakend))
                    continue;

                List<LinkedPair> lowerPairs = null;

                final SvVarData lowerSV = lowerBreakend.getSV();

                boolean lowerValidAP = mUseAlleleJCNs && mClusterJcnLimits.hasValidAlleleJcnData(
                        getClusterChrBreakendIndex(lowerBreakend), alleleJCNs);

                double lowerJcn = mLinkAllocator.getUnlinkedBreakendCount(lowerBreakend);

                int skippedNonAssembledIndex = -1; // the first index of a non-assembled breakend after the current one

                for (int j = i+1; j < breakendList.size(); ++j)
                {
                    final SvBreakend upperBreakend = breakendList.get(j);

                    if(belowJcnThreshold(upperBreakend.getSV()))
                        continue;

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

                    if((lowerSV.type() == INF || upperSV.type() == INF) && distance <= 1)
                    {
                        if(!jcnExceedsMajorAlleleJcn(lowerBreakend) || !jcnExceedsMajorAlleleJcn(upperBreakend))
                            continue;
                    }

                    LinkedPair newPair = new LinkedPair(lowerBreakend, upperBreakend);

                    // make note of adjacent facing breakends since these are prioritised in the adjacent pairs rule
                    boolean areAdjacent = false;
                    if(j == i + 1)
                    {
                        if(lowerBreakend.getDBLink() == null || lowerBreakend.getDBLink() != upperBreakend.getDBLink())
                            areAdjacent = true;
                    }
                    else if(j == i + 2)
                    {
                        // factor in DBs in between
                        SvBreakend middleBreakend = breakendList.get(i + 1);

                        if((lowerBreakend.getDBLink() != null && lowerBreakend.getDBLink() == middleBreakend.getDBLink())
                        || (upperBreakend.getDBLink() != null && upperBreakend.getDBLink() == middleBreakend.getDBLink()))
                            areAdjacent = true;
                    }

                    if(areAdjacent)
                    {
                        mAdjacentPairs.add(newPair);

                        if(copyNumbersEqual(lowerJcn, mLinkAllocator.getUnlinkedBreakendCount(upperBreakend)))
                            mAdjacentMatchingPairs.add(newPair);
                    }

                    if(lowerPairs == null)
                    {
                        lowerPairs = Lists.newArrayList();
                        mSvBreakendPossibleLinks.put(lowerBreakend, lowerPairs);
                    }

                    lowerPairs.add(newPair);

                    List<LinkedPair> upperPairs = mSvBreakendPossibleLinks.get(upperBreakend);

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

                    if(lowerValidAP && mClusterJcnLimits.hasValidAlleleJcnData(getClusterChrBreakendIndex(upperBreakend), alleleJCNs))
                    {
                        int breakendIndex = getClusterChrBreakendIndex(upperBreakend);
                        double clusterAP = alleleJCNs.get(breakendIndex).clusterJcn();

                        if(clusterAP < CLUSTER_ALLELE_JCN_MIN)
                        {
                            // this lower breakend cannot match with anything further upstream
                            LNX_LOGGER.trace("breakends lower({}: {}) limited at upper({}: {}) with clusterAP({})",
                                    i, lowerBreakend.toString(), j, upperBreakend.toString(), formatJcn(clusterAP));

                            break;
                        }
                    }
                }
            }
        }
    }

    private void checkIsComplexDupSV(SvBreakend lowerJcnBreakend, SvBreakend higherJcnBreakend)
    {
        // check if the lower JCN SV connects to both ends of another SV to replicate it
        final SvVarData var = lowerJcnBreakend.getSV();

        if(var.jcn() < 1) // maintain a minimum to avoid JCN comparison issues for lower values
            return;

        if(var.isSglBreakend() || var.type() == DEL || higherJcnBreakend.getSV().isSglBreakend())
            return;

        if(mComplexDupCandidates.keySet().contains(var))
            return;

        if(mLinkAllocator.getSvConnections().get(var) == null)
            return;

        final SvVarData otherSV = higherJcnBreakend.getSV();

        if(var.jcn() > otherSV.jcn() || var.jcnMin() * 2 > otherSV.jcnMax())
            return;

        if(!jcnMatchForSplits(var.jcn(), var.jcnUncertainty(), otherSV.jcn(), otherSV.jcnUncertainty()))
            return;

        // check whether the other breakend satisfies the same ploidy comparison criteria
        SvBreakend otherBreakend = var.getBreakend(!lowerJcnBreakend.usesStart());

        final List<SvBreakend> breakendList = mChrBreakendMap.get(otherBreakend.chromosome());

        boolean traverseUp = otherBreakend.orientation() == -1;
        int index = getClusterChrBreakendIndex(otherBreakend);

        while(true)
        {
            index += traverseUp ? 1 : -1;

            if(index < 0 || index >= breakendList.size())
                break;

            final SvBreakend breakend = breakendList.get(index);

            if(breakend == lowerJcnBreakend)
                break;

            if (breakend.isAssembledLink())
                continue;

            if (breakend.orientation() == otherBreakend.orientation())
                break;

            if(Math.abs(breakend.position() - otherBreakend.position()) < getMinTemplatedInsertionLength(breakend, otherBreakend))
                continue;

            SvVarData otherSV2 = breakend.getSV();

            if(!jcnMatchForSplits(var.jcn(), var.jcnUncertainty(), otherSV2.jcn(), otherSV2.jcnUncertainty()))
                return;

            List<LinkedPair> links = Lists.newArrayList(LinkedPair.from(lowerJcnBreakend, higherJcnBreakend),
                    LinkedPair.from(otherBreakend, breakend));

            if(otherSV == otherSV2)
            {
                logInfo(String.format("identified complex dup(%s %s) ploidy(%s) vs SV(%s) ploidy(%s)",
                        var.posId(), var.type(), formatJcn(var.jcn()), otherSV.id(), formatJcn(otherSV.jcn())));
            }
            else
            {
                logInfo(String.format("identified complex dup(%s %s) ploidy(%s) vs SVs(%s & %s) ploidy(%s & %s)",
                        var.posId(), var.type(), formatJcn(var.jcn()), otherSV.id(), otherSV2.id(),
                        formatJcn(otherSV.jcn()), formatJcn(otherSV2.jcn())));
            }

            mComplexDupCandidates.put(var, links);

            break;
        }
    }

    private boolean alreadyLinkedBreakend(final SvBreakend breakend)
    {
        // assembled links have already been added to chains prior to determining remaining possible links
        // so these need to be excluded unless their replication count allows them to be used again
        return breakend.isAssembledLink() && mLinkAllocator.getUnlinkedBreakendCount(breakend) == 0;
    }

    private void checkDoubleMinuteChains()
    {
        // if there is a single chain which contains all DM SVs, attempt to close the chain
        if(mDoubleMinuteSVs.isEmpty())
            return;

        // at the moment this logic is problematic and can divide chains repeatedly
        // reconcileChains(mChains, true, mLinkAllocator.getNextChainId(), true);

        // search for a chain which can be closed if it contains all the DM SVs from one of the closed DM chains (found by DM Finder)
        final List<SvChain> dmChains = mDoubleMinuteSVs.get(0).getCluster().getDoubleMinuteChains();
        final List<SvChain> unmatchedDmChains = Lists.newArrayList(dmChains);

        for(SvChain chain : mChains)
        {
            // matches DM chain
            boolean matchesDmChain = false;
            for(SvChain dmChain : dmChains)
            {
                if(chain.hasMatchingSVs(dmChain))
                {
                    matchesDmChain = true;
                    break;
                }
            }

            if(matchesDmChain && chain.couldCloseChain())
            {
                SvBreakend chainStart = chain.getOpenBreakend(true);
                SvBreakend chainEnd = chain.getOpenBreakend(false);

                if(chainEnd.getSV() == chainStart.getSV() && chainStart != chainEnd)
                {
                    chain.closeChain();
                }
                else
                {
                    chain.closeChain(LR_METHOD_DM_CLOSE, mLinkAllocator.getLinkIndex());
                    LNX_LOGGER.debug("cluster({}) closed DM chain({})", mClusterId, chain.id());
                }
            }

            for(SvChain dmChain : dmChains)
            {
                if(identicalChain(chain, dmChain, false, true))
                {
                    chain.setDoubleMinute(true);
                    unmatchedDmChains.remove(dmChain);
                    break;
                }
            }
        }

        if(!unmatchedDmChains.isEmpty())
        {
            LNX_LOGGER.debug("cluster({}) added pre-formed DM chain", mClusterId);
            unmatchedDmChains.forEach(x -> mChains.add(x));
        }
    }

    private void checkChains()
    {
        for(final SvChain chain : mChains)
        {
            if(!checkIsValid(chain))
            {
                LNX_LOGGER.error("cluster({}) has invalid chain({})", mClusterId, chain.id());
                chain.logLinks();
                mIsValid = false;
            }
        }

        // check no 2 chains have the same link reference
        for(int i = 0; i < mChains.size() - 1; ++i)
        {
            final SvChain chain1 = mChains.get(i);

            for(int j = i + 1; j < mChains.size(); ++j)
            {
                final SvChain chain2 = mChains.get(j);

                for(final LinkedPair pair : chain2.getLinkedPairs())
                {
                    if(chain1.getLinkedPairs().stream().anyMatch(x -> x == pair))
                    {
                        LNX_LOGGER.error("cluster({}) chain({}) and chain({}) share pair({})",
                                mClusterId, chain1.id(), chain2.id(), pair.toString());
                        mIsValid = false;
                    }
                }
            }
        }
    }

    protected void logInfo(final String message)
    {
        mDiagnostics.addMessage(message);
        LNX_LOGGER.debug(message);
    }

    public void setLogVerbose(boolean toggle)
    {
        mLogVerbose = toggle;
        setRunValidation(toggle);
    }

    private void enableLogVerbose()
    {
        if(!mLogVerbose)
            return;

        mLogLevel = LNX_LOGGER.getLevel();
        Configurator.setRootLevel(TRACE);
    }

    private void disableLogVerbose()
    {
        if(!mLogVerbose)
            return;

        // restore logging
        Configurator.setRootLevel(mLogLevel);
    }

}
