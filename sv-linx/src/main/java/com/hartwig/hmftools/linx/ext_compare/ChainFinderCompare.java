package com.hartwig.hmftools.linx.ext_compare;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxConfig.formOutputPath;
import static com.hartwig.hmftools.linx.ext_compare.CfSvChainData.CF_DATA_DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class ChainFinderCompare
{
    private final String mDataDirectory;
    private CfSampleData mSampleData;

    public static final String CHAIN_FINDER_SAMPLE_DATA_DIR = "chain_finder_data_dir";

    public ChainFinderCompare(final CommandLine cmd)
    {
        mDataDirectory = formOutputPath(cmd.getOptionValue(CHAIN_FINDER_SAMPLE_DATA_DIR));
        mSampleData = null;
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(CHAIN_FINDER_SAMPLE_DATA_DIR, true, "Directory containing each sample's ChainFinder run results");
    }

    public void processSample(final String sampleId, final List<SvVarData> svDataList, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        final String resultsDir = mDataDirectory + sampleId + File.separator + "results" + File.separator + sampleId + File.separator;

        if(!Files.exists(Paths.get(resultsDir)))
        {
            LNX_LOGGER.error("sample({}) invalid file path({}) to results data", sampleId, resultsDir);
            return;
        }

        final String chainDataFile =  sampleId + "_chains_final.txt";

        if(!Files.exists(Paths.get(resultsDir + chainDataFile)))
        {
            LNX_LOGGER.error("sample({}) invalid file path to chain results data({} + {})", sampleId, resultsDir, chainDataFile);
            return;
        }

        mSampleData = new CfSampleData(sampleId);

        final List<CfBreakendData> breakendDataList = loadChainFinderData(sampleId, resultsDir + chainDataFile);

        mSampleData.UnchainedSvList.addAll(svDataList.stream()
                .filter(x -> !x.isSglBreakend()) // SGLs aren't handled by CF and so are considered out of scope
                .collect(Collectors.toList()));

        mapChainSvData(breakendDataList, chrBreakendMap);

        reportSVChaining();

        reportDeletionBridges();
    }

    private void mapChainSvData(final List<CfBreakendData> breakendDataList, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
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
            CfChain chain = mSampleData.Chains.get(cfSvData.ChainId);

            if(chain == null)
            {
                chain = new CfChain(cfSvData.ChainId);
                mSampleData.Chains.put(cfSvData.ChainId, chain);
            }

            if(!mapToLinx(cfSvData, chrBreakendMap))
            {
                LNX_LOGGER.warn("CF sv({}) not matched with actual SV", cfSvData);
                continue;
            }

            mSampleData.UnchainedSvList.remove(cfSvData.getSvData());
        }
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

    private static final int CF_MIN_CHAIN_COUNT = 3;

    private void reportSVChaining()
    {
        int linxChainedOnly = (int)mSampleData.UnchainedSvList.stream()
                .filter(x -> x.getCluster().getSvCount() >= CF_MIN_CHAIN_COUNT)
                .count();

        int cfChainedOnly = (int)mSampleData.CfSvList.stream().filter(x -> x.getSvData().getCluster().getSvCount() == 1).count();

        LNX_LOGGER.info("sample({}) SVs chaining matched({}) cfOnly({}) linxOnly({})",
                mSampleData.SampleId, mSampleData.CfSvList.size(), cfChainedOnly, linxChainedOnly);
    }

    public void reportClusteringStats()
    {

    }

    public void reportChainingStats()
    {

    }

    public void reportDeletionBridges()
    {


    }

    private void writeSvMatchData(final SvVarData var, final CfSvChainData cfSvData)
    {


    }

    private final List<CfBreakendData> loadChainFinderData(final String sampleId, final String chainDataFile)
    {
        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(chainDataFile));

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
