package com.hartwig.hmftools.linx.ext_compare;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxOutput.ITEM_DELIM;
import static com.hartwig.hmftools.linx.ext_compare.CfDbMatchType.LINX_ONLY;
import static com.hartwig.hmftools.linx.ext_compare.CfSvData.CF_DATA_DELIMITER;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.asStr;
import static com.hartwig.hmftools.linx.types.LinxConstants.NO_DB_MARKER;
import static com.hartwig.hmftools.linx.types.SvVarData.CR_DELIM;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.types.DbPair;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.Nullable;

public class ChainFinderCompare
{
    private final String mDataDirectory;
    private final boolean mUseSampleDirectories;
    private final boolean mLogDbData;

    private CfSampleData mSampleData;
    private BufferedWriter mSvDataWriter;
    private BufferedWriter mClusterChainOverlapWriter;

    public static final String CHAIN_FINDER_DATA_DIR = "chain_finder_data_dir";
    private static final String USE_SAMPLE_DIRECTORIES = "use_sample_dir";

    // chain-finder won't form chains of only 2 SVs
    private static final int CF_MIN_CHAIN_COUNT = 3;

    public ChainFinderCompare(final String outputDir, final CommandLine cmd)
    {
        mDataDirectory = checkAddDirSeparator(cmd.getOptionValue(CHAIN_FINDER_DATA_DIR));
        mUseSampleDirectories = cmd.hasOption(USE_SAMPLE_DIRECTORIES);
        mSampleData = null;
        mSvDataWriter = null;
        mClusterChainOverlapWriter = null;
        mLogDbData = true;

        initialiseSvDataWriter(outputDir);
        initialiseClusterChainOverlapWriter(outputDir);
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(CHAIN_FINDER_DATA_DIR, true, "Root directory for ChainFinder sample run results");
        options.addOption(USE_SAMPLE_DIRECTORIES, false, "All sample ChainFinder results in a single directory");
    }

    public void close()
    {
        closeBufferedWriter(mSvDataWriter);
        closeBufferedWriter(mClusterChainOverlapWriter);
    }

    public void processSample(
            final String sampleId, final List<SvVarData> svDataList,
            final List<SvCluster> clusters, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        final String resultsDir = mUseSampleDirectories ?
                mDataDirectory + sampleId + File.separator + "results" + File.separator + sampleId + File.separator : mDataDirectory;

        if(!Files.exists(Paths.get(resultsDir)))
        {
            LNX_LOGGER.error("sample({}) invalid path({}) to chain-finder results data", sampleId, resultsDir);
            return;
        }

        // check whether any chaining results have been written for this sample - it's possible for CF to not find any chains
        final String chainDataFile =  sampleId + "_chains_final.txt";

        mSampleData = new CfSampleData(sampleId);

        // set clustered SV distances
        // final List<SvVarData> unclusteredSVs = clusters.stream()
        //        .filter(x -> x.getSvCount() == 1).map(x -> x.getSV(0)).collect(Collectors.toList());

        clusters.forEach(x -> mSampleData.setSvClusterDistances(x));

        if(!Files.exists(Paths.get(resultsDir + chainDataFile)))
        {
            LNX_LOGGER.info("sample({}) invalid chain-finder results file({} + {})", sampleId, resultsDir, chainDataFile);
            writeUnchainedSVs(svDataList);
            return;
        }

        final List<CfBreakendData> breakendDataList = loadChainFinderData(sampleId, resultsDir + chainDataFile);

        if(breakendDataList.isEmpty())
            return;

        mSampleData.UnchainedSvList.addAll(svDataList.stream()
                .filter(x -> !x.isSglBreakend()) // SGLs aren't handled by CF and so are considered out of scope
                .collect(Collectors.toList()));

        mapChainSvData(breakendDataList, chrBreakendMap);

        mSampleData.setUnclusteredDistances(chrBreakendMap);

        for(CfSvData cfSvData : mSampleData.CfSvList)
        {
            if(cfSvData.getSvData() == null)
                continue;

            writeMatchData(cfSvData.getSvData(), cfSvData);
        }

        // write out any SV clustered by Linx if it's not chained by CF
        writeUnchainedSVs(mSampleData.UnchainedSvList);

        mSampleData.ClusterChainOverlaps.forEach(x -> writeClusterChainOverlapData(x));
    }

    private void mapChainSvData(final List<CfBreakendData> breakendDataList, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        // first correct the NaN chromosomes which should be either X or Y
        for(CfBreakendData breakend : breakendDataList)
        {
            if(breakend.Chromosome.equals("NaN"))
            {
                List<SvBreakend> breakends = chrBreakendMap.get("X");
                if(breakends != null && breakends.stream().anyMatch(x -> x.position() == breakend.Position))
                {
                    breakend.Chromosome = "X";
                }

                breakends = chrBreakendMap.get("Y");
                if(breakends != null && breakends.stream().anyMatch(x -> x.position() == breakend.Position))
                {
                    breakend.Chromosome = "Y";
                }
            }
        }

        // put chained breakend data pairs together
        int index = 0;
        while(index < breakendDataList.size())
        {
            final CfBreakendData firstBreakend = breakendDataList.get(index);

            boolean found = false;
            for(int j = index + 1; j < breakendDataList.size(); ++j)
            {
                if(breakendDataList.get(j).RearrangeId == firstBreakend.RearrangeId)
                {
                    mSampleData.CfSvList.add(new CfSvData(new CfBreakendData[] { firstBreakend, breakendDataList.get(j)} ));
                    breakendDataList.remove(j);
                    found = true;
                    break;
                }
            }

            if(found)
                breakendDataList.remove(index);
            else
                ++index;
        }

        if(!breakendDataList.isEmpty())
        {
            LNX_LOGGER.warn("sampleId({}) has {} unmatched breakends", mSampleData.SampleId, breakendDataList.size());
        }

        // find and link to corresponding SVs
        for(CfSvData cfSvData : mSampleData.CfSvList)
        {
            if(!mapToLinx(cfSvData, chrBreakendMap))
            {
                // likely due to Linx filtering
                LNX_LOGGER.debug("CF sv({}) not matched with actual SV", cfSvData);
                continue;
            }

            mSampleData.processNewSV(cfSvData);
        }

        mSampleData.setDeletionBridgeData();
        mSampleData.setChainLinks();
    }

    private boolean mapToLinx(final CfSvData cfSvData, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        final List<SvBreakend> breakendList = chrBreakendMap.get(cfSvData.Chromosomes[SE_START]);

        if(breakendList != null)
        {
            final SvBreakend matchingBreakend = breakendList.stream().filter(x -> cfSvData.matches(x.getSV())).findFirst().orElse(null);

            if(matchingBreakend != null)
            {
                cfSvData.setSvData(matchingBreakend.getSV());
                return true;
            }
        }

        return false;
    }

    private void reportSVChaining()
    {
        int linxChainedOnly = (int)mSampleData.UnchainedSvList.stream()
                .filter(x -> x.getCluster().getSvCount() >= CF_MIN_CHAIN_COUNT)
                .count();

        int cfChainedOnly = (int)mSampleData.CfSvList.stream().filter(x -> x.getSvData().getCluster().getSvCount() == 1).count();

        LNX_LOGGER.info("sample({}) SVs chaining matched({}) cfOnly({}) linxOnly({})",
                mSampleData.SampleId, mSampleData.CfSvList.size(), cfChainedOnly, linxChainedOnly);
    }

    private void writeUnchainedSVs(final List<SvVarData> svList)
    {
        svList.stream()
                .filter(x -> x.getCluster().getSvCount() - x.getCluster().getSglBreakendCount() >= 2)
                .forEach(x -> writeMatchData(x, null));
    }

    private void initialiseSvDataWriter(final String outputDir)
    {
        try
        {
            final String outputFileName = outputDir + "LNX_CHAIN_FINDER_SVS.csv";
            mSvDataWriter = createBufferedWriter(outputFileName, false);

            mSvDataWriter.write("SampleId,SvId,ClusterId,ClusterCount,SglCount,ResolvedType,Type,ClusterReason,ProxDistance");
            mSvDataWriter.write(",CfSvId,CfChainId,CfChainCount,SharedSvCount,OverlapGroupId,CfProxDistance");
            mSvDataWriter.write(",DbMatchedStart,DbMatchedEnd,CfLinkSvIdStart,CfLinkSvIdEnd,CfLinkLenStart,CfLinkLenEnd");

            // written again for convenience
            mSvDataWriter.write(",ChrStart,PosStart,OrientStart,ArmStart,ChrEnd,PosEnd,OrientEnd,ArmEnd");
            mSvDataWriter.write(",ClusteringDetails,Ploidy,ChainId,ChainCount");

            if(mLogDbData)
            {
                mSvDataWriter.write(",LnxDbLenStart,LnxDbLenEnd,CfDbLenStart,CfDbLenEnd");
            }

            mSvDataWriter.newLine();
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("error writing chain-finder SV data file: {}", e.toString());
        }
    }

    private String parseMainClusterReason(final String clusterReason)
    {
        if(clusterReason.isEmpty())
            return "NONE";

        final String[] reasons = clusterReason.split(ITEM_DELIM, -1);
        final String firstReason = reasons[0].split(CR_DELIM)[0];

        return firstReason;
    }

    private void writeMatchData(final SvVarData var, @Nullable final CfSvData cfSvData)
    {
        try
        {
            Integer proximityDistance = mSampleData.SvProximityDistance.get(var);

            mSvDataWriter.write(String.format("%s,%d,%d,%d,%d,%s,%s,%s,%d",
                mSampleData.SampleId, var.id(), var.getCluster().id(), var.getCluster().getSvCount(),
                    var.getCluster().getSglBreakendCount(), var.getCluster().getResolvedType(), var.type(),
                    parseMainClusterReason(var.getClusterReason()), proximityDistance != null ? proximityDistance : -1));

            final DbPair[] dbLinks = new DbPair[] { var.getDBLink(true), var.getDBLink(false) };

            if(cfSvData != null)
            {
                final CfChain chain = mSampleData.Chains.get(cfSvData.ChainId);
                final CfChainClusterOverlap clusterOverlap = mSampleData.findClusterChainOverlap(chain);

                if(chain == null || clusterOverlap == null)
                {
                    LNX_LOGGER.error("cfSvid({}) invalid chain or overlap group", cfSvData);
                    return;
                }

                final CfDbMatchType[] dbMatchTypes = cfSvData.getDeletionBridgeMatchTypes();

                mSvDataWriter.write(String.format(",%d,%d,%d,%d,%d,%d,%s,%s",
                        cfSvData.RearrangeId, cfSvData.ChainId, chain.SVs.size(), chain.getSharedSvCount(cfSvData),
                        clusterOverlap.Id, chain.findMinBreakendDistance(cfSvData), dbMatchTypes[SE_START], dbMatchTypes[SE_END]));

                final CfLink[] cfLinks = cfSvData.getChainLinks();

                mSvDataWriter.write(String.format(",%d,%d,%d,%d",
                        cfLinks[SE_START] != null ? cfLinks[SE_START].getOtherSv(cfSvData).getSvData().id() : -1,
                        cfLinks[SE_END] != null ? cfLinks[SE_END].getOtherSv(cfSvData).getSvData().id() : -1,
                        cfLinks[SE_START] != null ? cfLinks[SE_START].length() : -1,
                        cfLinks[SE_END] != null ? cfLinks[SE_END].length() : -1));
            }
            else
            {
                mSvDataWriter.write(String.format(",-1,-1,0,0,-1,-1,%s,%s,-1,-1,-1,-1",
                        dbLinks[SE_START] != null ? LINX_ONLY : CfDbMatchType.NONE,
                        dbLinks[SE_END] != null ? LINX_ONLY : CfDbMatchType.NONE));

            }

            mSvDataWriter.write(String.format(",%s,%d,%d,%s,%s,%d,%d,%s",
                    var.chromosome(true), var.position(true), var.orientation(true), asStr(var.arm(true)),
                    var.chromosome(false), var.position(false), var.orientation(false), asStr(var.arm(false))));

            final SvChain chain = var.getCluster().findChain(var);

            mSvDataWriter.write(String.format(",%s,%.2f,%d,%d",
                    var.getClusterReason(), var.jcn(),
                    chain != null ? chain.id() : -1, chain != null ? chain.getSvCount() : 0));

            if(mLogDbData)
            {
                mSvDataWriter.write(String.format(",%d,%d,%d,%d",
                        dbLinks[SE_START] != null ? dbLinks[SE_START].length() : NO_DB_MARKER,
                        dbLinks[SE_END] != null ? dbLinks[SE_END].length() : NO_DB_MARKER,
                        cfSvData != null ? cfSvData.getDbLengths()[SE_START] : NO_DB_MARKER,
                        cfSvData != null ? cfSvData.getDbLengths()[SE_END] : NO_DB_MARKER));
            }

            mSvDataWriter.newLine();
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("error writing chain-finder SV data file: {}", e.toString());
        }
    }

    private void initialiseClusterChainOverlapWriter(final String outputDir)
    {
        try
        {
            final String outputFileName = outputDir + "LNX_CHAIN_FINDER_GROUPS.csv";
            mClusterChainOverlapWriter = createBufferedWriter(outputFileName, false);

            mClusterChainOverlapWriter.write(CfChainClusterOverlap.header());
            mClusterChainOverlapWriter.newLine();
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("error writing chain-finder cluster-chain overlap data file: {}", e.toString());
        }
    }

    private void writeClusterChainOverlapData(final CfChainClusterOverlap overlapData)
    {
        try
        {
            mClusterChainOverlapWriter.write(String.format("%s,%s", mSampleData.SampleId, overlapData.toCsv()));
            mClusterChainOverlapWriter.newLine();
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("error writing chain-finder cluster-chain overlap data file: {}", e.toString());
        }
    }

    private List<CfBreakendData> loadChainFinderData(final String sampleId, final String chainDataFile)
    {
        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(chainDataFile));

            if(lines.size() <= 2)
                return Lists.newArrayList();

            final Map<String,Integer> fielIndexMap = createFieldsIndexMap(lines.get(0), CF_DATA_DELIMITER);
            lines.remove(0); // header
            lines.remove(lines.size() - 1); // summary row

            return lines.stream().map(x -> CfBreakendData.fromData(x, fielIndexMap)).collect(Collectors.toList());
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to load chain data file({}): {}", chainDataFile.toString(), e.toString());
            return null;
        }

    }
}
