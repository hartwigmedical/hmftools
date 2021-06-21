package com.hartwig.hmftools.linx.chaining;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

public class ChainDiagnostics
{
    private List<String> mLogMessages;
    private String mOutputDir;
    private int mMaxClusterSize;
    private BufferedWriter mFileWriter;

    private final List<SvVarData> mInitialComplexDup;
    private final List<SvVarData> mInitialFoldbacks;
    private int mClusterCount;
    private int mUnlinkedBreakendCount;
    private int mUnlinkedSvCount;
    private int mWarnings;
    private int mLastProgressIndex;

    // references from chain finder
    private int mClusterId;
    private String mSampleId;
    private boolean mHasReplication;
    private final List<SvChain> mChains;
    private final List<SvChain> mUniqueChains;
    private final SvChainConnections mSvConnections;
    private final List<ChainState> mSvCompletedConnections;
    private final Map<SvBreakend, List<LinkedPair>> mSvBreakendPossibleLinks;
    private final List<SvVarData> mDoubleMinuteSVs;
    private final List<LinkedPair> mUniquePairs;

    public ChainDiagnostics(final SvChainConnections svConnMap, final List<ChainState> svCompleteConns,
            final List<SvChain> chains, final List<SvChain> uniqueChains,
            final Map<SvBreakend, List<LinkedPair>> svBreakendPossibleLinks,
            final List<SvVarData> doubleMinuteSVs, final List<LinkedPair> uniquePairs)
    {
        mLogMessages = Lists.newArrayList();
        mInitialComplexDup = Lists.newArrayList();
        mInitialFoldbacks = Lists.newArrayList();
        mOutputDir = null;
        mFileWriter = null;
        mMaxClusterSize = 0;
        mUnlinkedSvCount = 0;
        mUnlinkedBreakendCount = 0;
        mWarnings = 0;
        mClusterCount = 0;
        mLastProgressIndex = 0;
        mSampleId = "";

        mSvConnections = svConnMap;
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

    public void setSampleId(final String sampleId) { mSampleId = sampleId; }

    public int unlinkedSvCount() { return mUnlinkedSvCount; }

    public void clear()
    {
        mLogMessages.clear();
        mWarnings = 0;
        mClusterCount = 0;
        mLastProgressIndex = 0;
    }

    public void initialise(int clusterId, boolean hasReplication)
    {
        clear();

        mClusterId = clusterId;
        mHasReplication = hasReplication;

        mClusterCount = mSvConnections.size();
    }

    public void setPriorityData(List<SvVarData> complexDups, final List<SvVarData> foldbacks)
    {
        mInitialFoldbacks.clear();
        mInitialFoldbacks.addAll(foldbacks);
        mInitialComplexDup.clear();
        mInitialComplexDup.addAll(complexDups);
    }

    public void addMessage(final String msg)
    {
        mLogMessages.add(String.format("INFO: %s", msg));
    }

    public void diagnoseChains()
    {
        if(mChains.isEmpty() || mClusterCount < 3 || mClusterCount > mMaxClusterSize)
            return;

        if(!mDoubleMinuteSVs.isEmpty()) // for now skip these
            return;

        List<ChainState> svConnections = Lists.newArrayList();
        svConnections.addAll(mSvCompletedConnections);
        svConnections.addAll(mSvConnections.values());

        int invalidBreakends = findMultiConnectionBreakends(svConnections);

        writeResults(svConnections, invalidBreakends);

        if(!LNX_LOGGER.isDebugEnabled())
            return;

        if(!mLogMessages.isEmpty())
        {
            LNX_LOGGER.debug("diagnostic log:");
            mLogMessages.stream().forEach(x -> LNX_LOGGER.debug(x));
        }
    }

    private int findMultiConnectionBreakends(List<ChainState> svConnections)
    {
        int invalidCount = 0;

        for(final ChainState svConn : svConnections)
        {
            if(svConn.uniqueConnections(true) <= 1 && svConn.uniqueConnections(false) <= 1)
                continue;

            // are multiple connections explained by foldbacks, DMs or complex DUPs?
            // ignore foldbacks formed by a single breakend and ignore DM SVs
            if(svConn.SV.isSingleBreakendFoldback() || mDoubleMinuteSVs.contains(svConn.SV))
                continue;

            // just take the end with the most connections or either if the same
            boolean useStart = svConn.uniqueConnections(true) > svConn.uniqueConnections(false);

            int connectionCount = svConn.uniqueConnections(useStart);

            int assembledLinks = svConn.SV.getAssembledLinkedPairs(useStart).size();

            if(connectionCount == assembledLinks)
                continue;

            int foldbackCons = 0;
            int compDupCons = 0;
            String otherSVs = "";

            for(final SvBreakend otherBreakend : svConn.getConnections(useStart))
            {
                final SvVarData otherVar = otherBreakend.getSV();

                if(mInitialFoldbacks.contains(otherVar))
                    ++foldbackCons;
                else if(mInitialComplexDup.contains(otherVar))
                    ++compDupCons;

                otherSVs = appendStr(otherSVs, otherVar.idStr(), ';');
            }

            if(connectionCount - foldbackCons - compDupCons > 1)
            {
                LNX_LOGGER.debug("cluster({} count={}) breakend({}) rep({}) conns({}) foldbacks({}) compDups({}) asmb({})",
                        mClusterId, mSvConnections.size(), svConn.SV.getBreakend(useStart).toString(),
                        svConn.Jcn, connectionCount, foldbackCons, compDupCons, assembledLinks);

                logCsv("MULTI_CONN", svConn.SV,
                        String.format("conns(%d) jcn(%s-%s-%s) foldbacks(%d) compDups(%d) asmb(%d) otherSVs(%s)",
                                connectionCount, formatJcn(svConn.MinJcn), formatJcn(svConn.Jcn), formatJcn(svConn.MaxJcn),
                                foldbackCons, compDupCons, assembledLinks, otherSVs));

                ++invalidCount;
            }
        }

        return invalidCount;
    }

    public void logCsv(final String type, final SvVarData var, final String otherInfo)
    {
        if(mMaxClusterSize == 0)
            return;

        LNX_LOGGER.info("CHAIN_DIAG: {},{},{},{},{}", type, mSampleId, mClusterId, var.id(), otherInfo);
    }

    public void chainingComplete()
    {
        mUnlinkedBreakendCount = mSvConnections.values().stream()
                .mapToInt(x -> (x.breakendCount(true) == 0 ? 1 : 0) + (x.breakendCount(false) == 0 ? 1 : 0)).sum();

        mUnlinkedSvCount = (int) mSvConnections.values().stream()
                .filter(x -> x.breakendCount(true) == 0 && x.breakendCount(false) == 0).count();

        if(!LNX_LOGGER.isDebugEnabled())
            return;

        if(mHasReplication)
        {
            double unlinkedJcnBreakends = mSvConnections.values().stream()
                .mapToDouble(x -> x.unlinked(true) + x.unlinked(false)).sum();

            double unlinkedJcnSVs = mSvConnections.values().stream().mapToDouble(x -> x.unlinked()).sum();

            LNX_LOGGER.debug("cluster({}) chaining finished: chains({} unique={} links={}) SVs({}) unlinked SVs({} jcn={}) breakends({} jcn={})",
                    mClusterId, mChains.size(), mUniqueChains.size(), mUniquePairs.size(), mClusterCount,
                    mUnlinkedSvCount, formatJcn(unlinkedJcnSVs), mUnlinkedBreakendCount, formatJcn(unlinkedJcnBreakends));
        }
        else
        {
            LNX_LOGGER.debug("cluster({}) chaining finished: chains({} links={}) SVs({}) unlinked SVs({}) breakends({})",
                    mClusterId, mUniqueChains.size(), mUniquePairs.size(), mClusterCount, mUnlinkedSvCount, mUnlinkedBreakendCount);
        }
    }

    private void writeResults(final List<ChainState> svConnections, int invalidBreakends)
    {
        if(mOutputDir == null)
            return;

        try
        {
            if (mFileWriter == null)
            {
                String outputFileName = mOutputDir;

                outputFileName += "LNX_CHAINS.csv";

                mFileWriter = createBufferedWriter(outputFileName, false);

                mFileWriter.write("SampleId,ClusterId,Replication,SvCount,JcnTotal,Chains,RepeatedChains,SGLs,Warnings");
                mFileWriter.write(",MaxJcn,UnlinksSVs,UnlinkedBEs,InvalidBEs,Foldbacks,CompDups");
                mFileWriter.newLine();
            }

            int sglCount = (int) svConnections.stream().filter(x -> x.SV.isSglBreakend()).count();

            double jcnTotal = svConnections.stream().mapToDouble(x -> x.Jcn).sum();

            double maxJcn = mHasReplication ? svConnections.stream().mapToDouble(x -> x.Jcn).max().getAsDouble() : 1;

            mFileWriter.write(String.format("%s,%d,%s,%d,%.1f,%d,%d,%d,%d",
                    mSampleId, mClusterId, mHasReplication, mClusterCount, jcnTotal,
                    mUniqueChains.size(), mChains.size() - mUniqueChains.size(), sglCount, mWarnings));

            mFileWriter.write(String.format(",%.1f,%d,%d,%d,%d,%d",
                    maxJcn, mUnlinkedSvCount, mUnlinkedBreakendCount, invalidBreakends,
                    mInitialFoldbacks.size(), mInitialComplexDup.size()));

            mFileWriter.newLine();
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing DM data: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mFileWriter);
    }

    public void checkProgress(int linkIndex)
    {
        if(!LNX_LOGGER.isDebugEnabled())
            return;

        if(!mHasReplication || mSvConnections.size() < 100)
            return;

        if(linkIndex >= mLastProgressIndex + 100)
        {
            mLastProgressIndex = linkIndex;
            LNX_LOGGER.debug("cluster({}) chaining progress: index({}) SVs(incomplete={} complete={}) partialChains({})",
                    mClusterId, linkIndex, mSvConnections.size(), mSvCompletedConnections.size(), mChains.size());
        }
    }

    public boolean checkHasValidState(int linkIndex)
    {
        // first check that the remaining possible links are supported by unlinked breakends
        boolean isValid = true;

        for(Map.Entry<SvBreakend,List<LinkedPair>> entry : mSvBreakendPossibleLinks.entrySet())
        {
            SvBreakend breakend = entry.getKey();

            List<SvBreakend> breakendList = Lists.newArrayList(); // mUnlinkedBreakendMap.get(breakend);

            if(breakendList == null)
            {
                LNX_LOGGER.error("cluster({}) runIndex({}): breakend({}) has {} possible pairs but no available breakends",
                        mClusterId, linkIndex, breakend.toString(), entry.getValue().size());

                isValid = false;
            }
        }

        return isValid;
    }

}
