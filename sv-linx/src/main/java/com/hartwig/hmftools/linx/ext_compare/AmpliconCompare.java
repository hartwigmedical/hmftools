package com.hartwig.hmftools.linx.ext_compare;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.checkOutputDir;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxOutput.SUBSET_DELIM;
import static com.hartwig.hmftools.linx.ext_compare.AmpliconClusterMatch.findAmplifyingClusters;
import static com.hartwig.hmftools.linx.ext_compare.AmpliconData.DATA_DELIM;
import static com.hartwig.hmftools.linx.ext_compare.AmpliconData.HIGH_SV_JCN_THRESHOLD;
import static com.hartwig.hmftools.linx.ext_compare.AmpliconData.extractSampleId;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNOT_DM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.sv.SvRegion;
import com.hartwig.hmftools.linx.drivers.DriverAmpData;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;

public class AmpliconCompare
{
    private final String mAmpliconDataFile;
    private String mSampleId;
    
    private final Map<String,List<AmpliconData>> mSampleAmpData;

    private BufferedWriter mResultsWriter;

    public static final String AMPLICON_DATA_FILE = "amplicon_data_file";

    public AmpliconCompare(final String outputDir, final CommandLine cmd)
    {
        mAmpliconDataFile = checkOutputDir(cmd.getOptionValue(AMPLICON_DATA_FILE));
        mSampleAmpData = Maps.newHashMap();
        mResultsWriter = null;
        mSampleId = "";

        initialiseResultsWriter(outputDir);
        loadData();
    }

    public void processSample(
            final String sampleId, final List<SvCluster> clusters, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        final List<AmpliconData> ampDataList = mSampleAmpData.get(sampleId);

        if(ampDataList == null)
            return;

        mSampleId = sampleId;

        for(final AmpliconData ampData : ampDataList)
        {
            checkClusterMatches(ampData, chrBreakendMap);
        }
    }

    private void checkClusterMatches(final AmpliconData ampData, final Map<String,List<SvBreakend>> chrBreakendMap)
    {
        final Set<SvCluster> matchingClusters = Sets.newHashSet();

        for(final SvRegion ampRegion : ampData.Regions)
        {
            final List<SvCluster> regionClusters =  findClustersMatchingRegion(ampData, ampRegion, chrBreakendMap);

            regionClusters.forEach(x -> matchingClusters.add(x));
        }

        if(matchingClusters.isEmpty())
        {
            writeMatchData(ampData, null);
        }
        else
        {
            matchingClusters.forEach(x -> writeMatchData(ampData, x));
        }
    }

    private List<SvCluster> findClustersMatchingRegion(
            final AmpliconData ampData, final SvRegion ampRegion, final Map<String,List<SvBreakend>> chrBreakendMap)
    {
        final List<SvCluster> candidateClusters = Lists.newArrayList();

        final List<SvBreakend> breakendList = chrBreakendMap.get(ampRegion.Chromosome);
        if(breakendList == null)
        {
            LNX_LOGGER.warn("amp(%s) has unmatched region(%s)", ampData, ampRegion);
            return candidateClusters;
        }

        // possible matches could be:
        // a link which overlaps or is within the AMP region
        // an SV with high JCN within or facing the AMP region
        Map<SvCluster, DriverAmpData> clusterMapData = findAmplifyingClusters(ampRegion, breakendList);
        candidateClusters.addAll(clusterMapData.keySet());

        return candidateClusters;
    }

    private void initialiseResultsWriter(final String outputDir)
    {
        try
        {
            final String outputFileName = outputDir + "LNX_AMPLICON_COMPARE.csv";
            mResultsWriter = createBufferedWriter(outputFileName, false);

            mResultsWriter.write("SampleId,AmpClusterId,AmpType,AmpSVs,AmpHighSVs,AmpChromosomes");
            mResultsWriter.write(",AmpMaxCN,AmpMaxJCN,AmpFoldbackCN,AmpRegions,AmpWidth");
            mResultsWriter.write(",ClusterId,ClusterCount,ResolvedType,ClusterDesc,ClusterHasDM");
            mResultsWriter.write(",ClusterHighSVs,ClusterChromosomes,ClusterMaxCN,ClusterMaxJCN,ClusterFoldbackJCN");

            mResultsWriter.newLine();
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("error writing amplicon comparison file: {}", e.toString());
        }
    }

    private void writeMatchData(final AmpliconData ampData, final SvCluster cluster)
    {
        try
        {
            mResultsWriter.write(String.format("%s", mSampleId));

            mResultsWriter.write(String.format(",%d,%s,%d,%d,%s",
                    ampData.ClusterId, ampData.Type, ampData.Breakends, ampData.HighBreakends, ampData.chromosomes()));

            mResultsWriter.write(String.format(",%.1f,%.1f,%.1f,%d,%d",
                    ampData.MaxCN, ampData.MaxJCN, ampData.FoldbackCN,
                    ampData.Regions.size(), ampData.Regions.stream().mapToLong(x -> x.baseLength()).sum()));

            if(cluster != null)
            {
                int highJcnCount = 0;
                double maxJcn = 0;
                double maxCN = 0;
                double foldbackJcn = 0;
                final List<String> chromosomes = Lists.newArrayList();

                for(final SvVarData var : cluster.getSVs())
                {
                    maxCN = max(maxCN, max(var.copyNumber(true), var.copyNumber(false)));
                    maxJcn = max(maxJcn, var.jcn());

                    if(var.isFoldback())
                    {
                        foldbackJcn += (var.isChainedFoldback() ? 0.5 : 1.0) * var.jcn();
                    }

                    if(var.jcn() >= HIGH_SV_JCN_THRESHOLD)
                    {
                        ++highJcnCount;

                        if(!chromosomes.contains(var.chromosome(true)))
                            chromosomes.add(var.chromosome(true));

                        if(!var.isSglBreakend() && !chromosomes.contains(var.chromosome(false)))
                            chromosomes.add(var.chromosome(false));
                    }
                }

                String highJcnChromosomes = appendStrList(chromosomes, SUBSET_DELIM);

                mResultsWriter.write(String.format(",%d,%d,%s,%s,%s",
                        cluster.id(), cluster.getSvCount(), cluster.getResolvedType(), cluster.getDesc(),
                        cluster.hasAnnotation(CLUSTER_ANNOT_DM)));

                mResultsWriter.write(String.format(",%d,%s,%.1f,%.1f,%.1f",
                        highJcnCount, highJcnChromosomes, maxCN, maxJcn, foldbackJcn));
            }
            else
            {
                mResultsWriter.write(",-1,0,,,false");
                mResultsWriter.write(",0,,0,0,0");
            }

            mResultsWriter.newLine();
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("error writing amplicon match data: {}", e.toString());
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(AMPLICON_DATA_FILE, true, "Input file for Amplicon data");
    }

    public void close()
    {
        closeBufferedWriter(mResultsWriter);
    }

    private void loadData()
    {
        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(mAmpliconDataFile));
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DATA_DELIM);
            lines.remove(0); // remove header

            AmpliconData currentAmpData = null;
            int ampCount = 0;

            for(String line : lines)
            {
                final String[] items = line.split(DATA_DELIM, -1);

                final String sampleId = extractSampleId(fieldsIndexMap, items);

                List<AmpliconData> ampDataList = mSampleAmpData.get(sampleId);

                if(ampDataList == null)
                {
                    ampDataList = Lists.newArrayList();
                    mSampleAmpData.put(sampleId, ampDataList);
                    currentAmpData = null;
                }

                AmpliconData ampData = AmpliconData.from(fieldsIndexMap, items);

                if(currentAmpData == null || currentAmpData.ClusterId != ampData.ClusterId)
                {
                    currentAmpData = ampData;
                    ampDataList.add(currentAmpData);
                    ++ampCount;
                }
                else
                {
                    currentAmpData.Regions.add(ampData.Regions.get(0));
                }
            }

            LNX_LOGGER.info("loaded amplicon data for {} samples, {} amplicons", mSampleAmpData.size(), ampCount);
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to load amplicon data file({}): {}", mAmpliconDataFile, e.toString());
            return;
        }

    }
}
