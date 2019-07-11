package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvChain;
import com.hartwig.hmftools.linx.types.SvChainState;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ChainDiagnostics
{
    private List<String> mLogMessages;
    private String mOutputDir;
    private int mMaxClusterSize;
    private BufferedWriter mFileWriter;

    private final List<SvVarData> mInitialComplexDup;
    private final List<SvVarData> mInitialFoldbacks;
    private int mUnlinkedBreakendCount;
    private int mUnlinkedSvCount;
    private int mCompleteSvCount;
    private int mMaxRepCount;
    private int mWarnings;

    // references from chain finder
    private int mClusterId;
    private boolean mHasReplication;
    private final List<SvChain> mChains;
    private final List<SvChain> mUniqueChains;
    private final Map<SvVarData,SvChainState> mSvConnectionsMap;
    private final List<SvChainState> mSvCompletedConnections;
    private final Map<SvBreakend, List<SvLinkedPair>> mSvBreakendPossibleLinks;
    private final List<SvVarData> mDoubleMinuteSVs;
    private final List<SvLinkedPair> mUniquePairs;

    public static final int LOG_TYPE_ERROR = 0;
    public static final int LOG_TYPE_WARN = 1;
    public static final int LOG_TYPE_INFO = 2;
    public static final int LOG_TYPE_VERBOSE = 3;

    private static final Logger LOGGER = LogManager.getLogger(ChainDiagnostics.class);

    public ChainDiagnostics(final Map<SvVarData,SvChainState> svConnMap, final List<SvChainState> svCompleteConns,
            final List<SvChain> chains, final List<SvChain> uniqueChains,
            final Map<SvBreakend, List<SvLinkedPair>> svBreakendPossibleLinks,
            final List<SvVarData> doubleMinuteSVs, final List<SvLinkedPair> uniquePairs)
    {
        mLogMessages = Lists.newArrayList();
        mInitialComplexDup = Lists.newArrayList();
        mInitialFoldbacks = Lists.newArrayList();
        mOutputDir = null;
        mFileWriter = null;
        mMaxClusterSize = 0;
        mUnlinkedSvCount = 0;
        mUnlinkedBreakendCount = 0;
        mCompleteSvCount = 0;
        mWarnings = 0;

        mSvConnectionsMap = svConnMap;
        mSvCompletedConnections = svCompleteConns;
        mChains = chains;
        mUniqueChains = uniqueChains;
        mSvBreakendPossibleLinks = svBreakendPossibleLinks;
        mDoubleMinuteSVs = doubleMinuteSVs;
        mUniquePairs = uniquePairs;
    }

    public void setOutputDir(final String dir, int maxLogSize)
    {
        mOutputDir = dir;
        mMaxClusterSize = maxLogSize;
    }

    public void initialise(int clusterId, boolean hasReplication)
    {
        mLogMessages.clear();
        mWarnings = 0;
        mUnlinkedSvCount = 0;
        mUnlinkedBreakendCount = 0;
        mCompleteSvCount = 0;

        mClusterId = clusterId;
        mHasReplication = hasReplication;

        if(mHasReplication)
        {
            mMaxRepCount = mSvConnectionsMap.values().stream().mapToInt(x -> x.MaxPloidy).max().getAsInt();
        }
        else
        {
            mMaxRepCount = 0;
        }
    }

    public void setPriorityData(List<SvVarData> complexDups, final List<SvVarData> foldbacks)
    {
        mInitialFoldbacks.clear();
        mInitialFoldbacks.addAll(foldbacks);
        mInitialComplexDup.clear();
        mInitialComplexDup.addAll(complexDups);
    }

    private static String logTypeStr(int type)
    {
        switch(type)
        {
            case LOG_TYPE_ERROR: return "Error";
            case LOG_TYPE_WARN: return "Warn";
            case LOG_TYPE_INFO: return "Info";
            default: return "Verbose";
        }
    }

    public void addMessage(int level, final String msg)
    {
        if(level <= LOG_TYPE_WARN)
            ++mWarnings;

        mLogMessages.add(String.format("%s: %s", logTypeStr(level), msg));
    }

    public void diagnoseChains(final String sampleId)
    {
        if(mSvConnectionsMap.keySet().size() < 4 || mSvConnectionsMap.keySet().size() > mMaxClusterSize)
            return;

        if(!mDoubleMinuteSVs.isEmpty()) // for now skip these
            return;

        int invalidBreakends = findMultiConnectionBreakends(sampleId);

        writeResults(sampleId, invalidBreakends);

        if(!LOGGER.isDebugEnabled())
            return;

        if(!mLogMessages.isEmpty())
        {
            LOGGER.debug("diagnostic log:");
            mLogMessages.stream().forEach(x -> LOGGER.debug(x));
        }
    }

    private int findMultiConnectionBreakends(final String sampleId)
    {
        if(mUniqueChains.size() > 2)
            return 0;

        List<SvBreakend> reportedBreakends = Lists.newArrayList();

        int invalidCount = 0;

        for(int i = 0; i < mUniquePairs.size(); ++i)
        {
            SvLinkedPair pair = mUniquePairs.get(i);

            for (int be = SE_START; be <= SE_END; ++be)
            {
                final SvBreakend breakend = pair.getBreakend(isStart(be));

                if (reportedBreakends.contains(breakend))
                    continue;

                reportedBreakends.add(breakend);

                // ignore foldbacks formed by a single breakend and ignore DM SVs
                if(breakend.getSV().isSingleBreakendFoldback() || mDoubleMinuteSVs.contains(breakend.getSV()))
                    continue;

                SvVarData otherVar = pair.getOtherBreakend(breakend).getSV();

                int connectionCount = 1;
                int foldbackCons = 0;
                int compDupCons = 0;

                if(mInitialFoldbacks.contains(otherVar))
                    ++foldbackCons;
                else if(mInitialComplexDup.contains(otherVar))
                    ++compDupCons;

                for (int j = i + 1; j < mUniquePairs.size(); ++j)
                {
                    SvLinkedPair pair2 = mUniquePairs.get(j);

                    SvVarData otherVar2 = null;

                    if(pair2.getBreakend(true) == breakend)
                        otherVar2 = pair2.getBreakend(false).getSV();
                    else if(pair2.getBreakend(false) == breakend)
                        otherVar2 = pair2.getBreakend(true).getSV();
                    else
                        continue;

                    ++connectionCount;

                    if(mInitialFoldbacks.contains(otherVar2))
                        ++foldbackCons;
                    else if(mInitialComplexDup.contains(otherVar2))
                        ++compDupCons;
                }

                if(connectionCount <= 1)
                    continue;

                int assembledLinks = breakend.getSV().getAssembledLinkedPairs(breakend.usesStart()).size();

                if(connectionCount == assembledLinks)
                    continue;

                if(connectionCount - foldbackCons - compDupCons > 1)
                {
                    int repCount = max(breakend.getSV().getReplicatedCount(), 1);

                    LOGGER.debug("cluster({} count={}) breakend({}) rep({}) conns({}) foldbacks({}) compDups({}) asmb({})",
                            mClusterId, mSvConnectionsMap.keySet().size(), breakend.toString(),
                            repCount, connectionCount, foldbackCons, compDupCons, assembledLinks);

                    ++invalidCount;
                }
            }
        }

        return invalidCount;
    }

    public void chainingComplete()
    {
        mUnlinkedBreakendCount = mSvConnectionsMap.values().stream()
                .mapToInt(x -> (x.breakendCount(true) == 0 ? 1 : 0) + (x.breakendCount(false) == 0 ? 1 : 0)).sum();

        mUnlinkedSvCount = (int) mSvConnectionsMap.values().stream().filter(x -> x.curentCount() == 0).count();

        mCompleteSvCount = mSvCompletedConnections.size();

        if(mHasReplication)
        {
            int unlinkedPloidyBreakends = mSvConnectionsMap.values().stream()
                .mapToInt(x -> x.unlinked(true) + x.unlinked(false)).sum();

            int unlinkedPloidySVs = mSvConnectionsMap.values().stream().mapToInt(x -> x.unlinked()).sum();

            LOGGER.debug("cluster({}) chaining finished: chains({} unique={} links={}) unlinked SVs({} unique={}) breakends({} reps={})",
                    mClusterId, mChains.size(), mUniqueChains.size(), mUniquePairs.size(), mUnlinkedSvCount, unlinkedPloidySVs,
                    mUnlinkedBreakendCount, unlinkedPloidyBreakends);
        }
        else
        {
            LOGGER.debug("cluster({}) chaining finished: chains({} links={}) unlinked SVs({}) breakends({})",
                    mClusterId, mUniqueChains.size(), mUniquePairs.size(), mUnlinkedSvCount, mUnlinkedBreakendCount);
        }
    }

    private void writeResults(final String sampleId, int invalidBreakends)
    {
        if(mOutputDir == null)
            return;

        try
        {
            if (mFileWriter == null)
            {
                String outputFileName = mOutputDir;

                outputFileName += "SVA_CHAINS.csv";

                mFileWriter = createBufferedWriter(outputFileName, false);

                mFileWriter.write("SampleId,ClusterId,Replication,SvCount,RepSvCount,Chains,ReplicatedChains,SGLs,Warnings");
                mFileWriter.write(",MaxRep,UnlinksSVs,UnlinkedBEs,InvalidBEs,Foldbacks,CompDups");
                mFileWriter.newLine();
            }

            int sglCount = (int) mSvConnectionsMap.keySet().stream().filter(SvVarData::isNullBreakend).count();

            // FIXME:
            int mIntialReplicatedSvCount = 0;

            mFileWriter.write(String.format("%s,%d,%s,%d,%d,%d,%d,%d,%d",
                    sampleId, mClusterId, mHasReplication, mSvConnectionsMap.keySet().size(), mIntialReplicatedSvCount,
                    mUniqueChains.size(), mChains.size() - mUniqueChains.size(), sglCount, mWarnings));

            mFileWriter.write(String.format(",%d,%d,%d,%d,%d,%d",
                    mMaxRepCount, mUnlinkedSvCount, mUnlinkedBreakendCount, invalidBreakends,
                    mInitialFoldbacks.size(), mInitialComplexDup.size()));

            mFileWriter.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing DM data: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mFileWriter);
    }

    public void checkProgress(int linkIndex)
    {
        if(!LOGGER.isDebugEnabled())
            return;

        /*
        if(!mHasReplication || mUniqueSVs.size() < 100)
            return;

        if((linkIndex % 100) == 0)
        {
            LOGGER.debug("cluster({}) chaining progress: SVs({}) partialChains({}) unlinked(SVs={} breakends={}) replicatedSVs({})",
                    mClusterId, mUniqueSVs.size(), mChains.size(), mUnlinkedReplicatedSVs.size(),
                    mUnlinkedBreakendMap.size(), mSvReplicationMap.size());
        }
        */
    }

    public boolean checkHasValidState(int linkIndex)
    {
        // first check that the remaining possible links are supported by unlinked breakends
        boolean isValid = true;

        for(Map.Entry<SvBreakend,List<SvLinkedPair>> entry : mSvBreakendPossibleLinks.entrySet())
        {
            SvBreakend breakend = entry.getKey();

            List<SvBreakend> breakendList = Lists.newArrayList(); // mUnlinkedBreakendMap.get(breakend);

            if(breakendList == null)
            {
                LOGGER.error("cluster({}) runIndex({}): breakend({}) has {} possible pairs but no available breakends",
                        mClusterId, linkIndex, breakend.toString(), entry.getValue().size());

                isValid = false;
            }
        }

        return isValid;
    }

}
