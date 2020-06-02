package com.hartwig.hmftools.linx.ext_compare;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxConfig.formOutputPath;
import static com.hartwig.hmftools.linx.analysis.ClusteringState.CR_FOLDBACKS;
import static com.hartwig.hmftools.linx.analysis.ClusteringState.CR_HOM_LOSS;
import static com.hartwig.hmftools.linx.analysis.ClusteringState.CR_LOH;
import static com.hartwig.hmftools.linx.analysis.ClusteringState.CR_LONG_DEL_DUP;
import static com.hartwig.hmftools.linx.analysis.ClusteringState.CR_LONG_INV;
import static com.hartwig.hmftools.linx.analysis.ClusteringState.CR_PROXIMITY;
import static com.hartwig.hmftools.linx.ext_compare.CfDbMatchType.LINX_ONLY;
import static com.hartwig.hmftools.linx.ext_compare.CfSvChainData.CF_DATA_DELIMITER;
import static com.hartwig.hmftools.linx.types.ChromosomeArm.asStr;
import static com.hartwig.hmftools.linx.types.SvConstants.NO_DB_MARKER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.Nullable;

public class ChainFinderCompare
{
    private final String mDataDirectory;
    private CfSampleData mSampleData;
    private BufferedWriter mSvDataWriter;

    public static final String CHAIN_FINDER_SAMPLE_DATA_DIR = "chain_finder_data_dir";

    // chain-finder won't form chains of only 2 SVs
    private static final int CF_MIN_CHAIN_COUNT = 3;

    public ChainFinderCompare(final String outputDir, final CommandLine cmd)
    {
        mDataDirectory = formOutputPath(cmd.getOptionValue(CHAIN_FINDER_SAMPLE_DATA_DIR));
        mSampleData = null;
        mSvDataWriter = null;

        initialiseSvDataWriter(outputDir);
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(CHAIN_FINDER_SAMPLE_DATA_DIR, true, "Directory containing each sample's ChainFinder run results");
    }

    public void close() { closeBufferedWriter(mSvDataWriter); }

    public void processSample(final String sampleId, final List<SvVarData> svDataList, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        final String resultsDir = mDataDirectory + sampleId + File.separator + "results" + File.separator + sampleId + File.separator;

        if(!Files.exists(Paths.get(resultsDir)))
        {
            LNX_LOGGER.error("sample({}) invalid path({}) to chain-finder results data", sampleId, resultsDir);
            return;
        }

        // check whether any chaining results have been written for this sample - it's possible for CF to not find any chains
        final String chainDataFile =  sampleId + "_chains_final.txt";

        if(!Files.exists(Paths.get(resultsDir + chainDataFile)))
        {
            LNX_LOGGER.debug("sample({}) invalid chain-finder results file({} + {})", sampleId, resultsDir, chainDataFile);
            writeUnchainedSVs(svDataList);
            return;
        }

        final List<CfBreakendData> breakendDataList = loadChainFinderData(sampleId, resultsDir + chainDataFile);

        if(breakendDataList.isEmpty())
            return;

        mSampleData = new CfSampleData(sampleId);

        mSampleData.UnchainedSvList.addAll(svDataList.stream()
                .filter(x -> !x.isSglBreakend()) // SGLs aren't handled by CF and so are considered out of scope
                .collect(Collectors.toList()));

        mapChainSvData(breakendDataList, chrBreakendMap);

        for(CfSvChainData cfSvData : mSampleData.CfSvList)
        {
            if(cfSvData.getSvData() == null)
                continue;

            writeMatchData(cfSvData.getSvData(), cfSvData);
        }

        // write out any SV clustered by Linx if it's not chained by CF
        writeUnchainedSVs(mSampleData.UnchainedSvList);
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
                    mSampleData.CfSvList.add(new CfSvChainData(new CfBreakendData[] { firstBreakend, breakendDataList.get(j)} ));
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
        for(CfSvChainData cfSvData : mSampleData.CfSvList)
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
    }

    private boolean mapToLinx(final CfSvChainData cfSvData, final Map<String, List<SvBreakend>> chrBreakendMap)
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

            mSvDataWriter.write("SampleId,SvId,ClusterId,ClusterCount,SglCount,ResolvedType,Type,ClusterReason");
            mSvDataWriter.write(",CfSvId,CfChainId,CfChainCount,SharedSvCount,DbMatchedStart,DbMatchedEnd");

            // written again for convenience
            mSvDataWriter.write(",ChrStart,PosStart,OrientStart,ArmStart,ChrEnd,PosEnd,OrientEnd,ArmEnd");
            mSvDataWriter.write(",ClusteringDetails,Ploidy,ChainId,ChainCount");

            // skip the detailed DB output for now
            // mSvDataWriter.write(",DbMatchedStart,DbMatchedEnd,LnxDbLenStart,LnxDbLenEnd,CfDbLenStart,CfDbLenEnd");

            // mSvDataWriter.write(",Foldback");

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

        final String[] reasons = clusterReason.split(";", -1);
        final String firstReason = reasons[0].split("_")[0];

        return firstReason;
    }

    private void writeMatchData(final SvVarData var, @Nullable final CfSvChainData cfSvData)
    {
        try
        {
            mSvDataWriter.write(String.format("%s,%d,%d,%d,%d,%s,%s,%s",
                mSampleData.SampleId, var.id(), var.getCluster().id(), var.getCluster().getSvCount(),
                    var.getCluster().getSglBreakendCount(), var.getCluster().getResolvedType(),
                    var.type(), parseMainClusterReason(var.getClusterReason())));

            final SvLinkedPair[] dbLinks = new SvLinkedPair[] { var.getDBLink(true), var.getDBLink(false) };

            if(cfSvData != null)
            {
                final CfChain chain = mSampleData.Chains.get(cfSvData.ChainId);
                final CfChainClusterOverlap clusterOverlap = mSampleData.LinkedClusterChains.get(cfSvData.clusterChainId());
                final CfDbMatchType[] dbMatchTypes = cfSvData.getDeletionBridgeMatchTypes();

                mSvDataWriter.write(String.format(",%d,%d,%d,%d,%s,%s",
                        cfSvData.RearrangeId, cfSvData.ChainId, chain.ChainSVs.size(), clusterOverlap.SharedSVs.size(),
                        dbMatchTypes[SE_START], dbMatchTypes[SE_END]));
                /*
                mSvDataWriter.write(String.format(",%d,%d,%d,%d,%s,%s,%d,%d,%d,%d",
                        dbLinks[SE_START] != null ? dbLinks[SE_START].length() : NO_DB_MARKER,
                        dbLinks[SE_END] != null ? dbLinks[SE_END].length() : NO_DB_MARKER,
                        cfSvData.getDbLengths()[SE_START], cfSvData.getDbLengths()[SE_END]));
                */
            }
            else
            {
                mSvDataWriter.write(String.format(",-1,-1,0,0,%s,%s",
                        dbLinks[SE_START] != null ? LINX_ONLY : CfDbMatchType.NONE,
                        dbLinks[SE_END] != null ? LINX_ONLY : CfDbMatchType.NONE));

                /*
                mSvDataWriter.write(String.format(",-1,-1,0,0,%s,%s,%d,%d,%d,%d",
                        dbLinks[SE_START] != null ? LINX_ONLY : CfDbMatchType.NONE,
                        dbLinks[SE_END] != null ? LINX_ONLY : CfDbMatchType.NONE,
                        dbLinks[SE_START] != null ? dbLinks[SE_START].length() : NO_DB_MARKER,
                        dbLinks[SE_END] != null ? dbLinks[SE_END].length() : NO_DB_MARKER, NO_DB_MARKER, NO_DB_MARKER));
                */
            }

            mSvDataWriter.write(String.format(",%s,%d,%d,%s,%s,%d,%d,%s",
                    var.chromosome(true), var.position(true), var.orientation(true), asStr(var.arm(true)),
                    var.chromosome(false), var.position(false), var.orientation(false), asStr(var.arm(false))));

            final SvChain chain = var.getCluster().findChain(var);

            mSvDataWriter.write(String.format(",%s,%.2f,%d,%d",
                    var.getClusterReason(), var.ploidy(),
                    chain != null ? chain.id() : -1, chain != null ? chain.getSvCount() : 0));

            mSvDataWriter.newLine();

            // mSvDataWriter.write(String.format(",%s,%.2f,%s",var.isFoldback()));
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("error writing chain-finder SV data file: {}", e.toString());
        }
    }

    private final List<CfBreakendData> loadChainFinderData(final String sampleId, final String chainDataFile)
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
