package com.hartwig.hmftools.svtools.mult_biopsy;

import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svtools.mult_biopsy.MultiBiopsyData.MATCH_TYPE_PARTIAL;
import static com.hartwig.hmftools.svtools.mult_biopsy.MultiBiopsyData.MATCH_TYPE_PRIVATE;
import static com.hartwig.hmftools.svtools.mult_biopsy.MultiBiopsyData.MATCH_TYPE_SHARED;
import static com.hartwig.hmftools.svtools.mult_biopsy.MultiBiopsyData.areMatched;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class MultipleBiopsyAnalyser
{
    private Map<String, List<String>> mPatientSampleIdsMap;
    private Map<String, List<MultiBiopsyData>> mSampleSvData;

    private BufferedWriter mSvWriter;
    private BufferedWriter mMergeWriter;
    private BufferedWriter mClusterOverlapWriter;
    private String mOutputDir;

    private static final Logger LOGGER = LogManager.getLogger(MultipleBiopsyAnalyser.class);

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        String outputDir = parseOutputDir(cmd);

        MultipleBiopsyAnalyser mbAnalyser = new MultipleBiopsyAnalyser();

        if(!mbAnalyser.loadData(cmd, outputDir))
            return;

        mbAnalyser.runAnalysis();

        LOGGER.info("multiple-biopsy analysis complete");
    }

    private static Options createBasicOptions()
    {
        final Options options = new Options();
        options.addOption(PATIENT_SAMPLE_IDS_FILE, true, "File mapping PatientIds to SampleIds file");
        options.addOption(SVS_INPUT_FILE, true, "LINX SVs file");
        addOutputDir(options);
        addLoggingOptions(options);
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    public MultipleBiopsyAnalyser()
    {
        mPatientSampleIdsMap = Maps.newHashMap();
        mSampleSvData = Maps.newHashMap();
        mSvWriter = null;
        mMergeWriter = null;
        mClusterOverlapWriter = null;
        mOutputDir = "";
    }

    private static String PATIENT_SAMPLE_IDS_FILE = "patient_ids_file";
    private static String SVS_INPUT_FILE = "sv_data_file";

    public boolean loadData(final CommandLine cmd, final String outputDir)
    {
        if (!loadPatientSampleData(cmd.getOptionValue(PATIENT_SAMPLE_IDS_FILE)))
            return false;

        if (!loadSampleSVData(cmd.getOptionValue(SVS_INPUT_FILE)))
            return false;

        mOutputDir = outputDir;

        return true;
    }

    private boolean loadPatientSampleData(final String filename)
    {
        if (filename.isEmpty() || !Files.exists(Paths.get(filename)))
            return false;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            // skip field names
            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("Empty patient sample IDs file({})", filename);
                return false;
            }

            int count = 0;
            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(",");

                if (items.length != 2)
                    continue;

                String patientId = items[0];
                String sampleId = items[1];

                List<String> sampleIds = mPatientSampleIdsMap.get(patientId);

                if (sampleIds == null)
                {
                    sampleIds = Lists.newArrayList();
                    mPatientSampleIdsMap.put(patientId, sampleIds);
                }

                sampleIds.add(sampleId);
                ++count;
            }

            LOGGER.info("loaded {} patient-sample ID records", count);
        } catch (IOException e)
        {
            LOGGER.error("Failed to load patient sample IDs file({}): {}", filename, e.toString());
        }

        return true;
    }

    private boolean loadSampleSVData(final String filename)
    {
        if (filename.isEmpty() || !Files.exists(Paths.get(filename)))
            return false;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("Empty MB SV data file({})", filename);
                return false;
            }

            int count = 0;
            String currentSampleId = "";
            List<MultiBiopsyData> mbDataList = null;
            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(",");

                if (items.length < 13)
                    continue;

                String sampleId = items[0];

                if (!currentSampleId.equals(sampleId))
                {
                    mbDataList = Lists.newArrayList();
                    mSampleSvData.put(sampleId, mbDataList);
                    currentSampleId = sampleId;
                }

                MultiBiopsyData mbData = new MultiBiopsyData(items);

                mbDataList.add(mbData);
                ++count;
            }

            LOGGER.info("loaded {} SV data records", count);
        }
        catch (IOException e)
        {
            LOGGER.error("failed to load patient sample IDs file({}): {}", filename, e.toString());
        }

        return true;
    }

    public void runAnalysis()
    {
        mapSamples();

        LOGGER.info("run complete");

        closeOutputFiles();
    }

    private void mapSamples()
    {
        for (Map.Entry<String, List<String>> entry : mPatientSampleIdsMap.entrySet())
        {
            String patientId = entry.getKey();
            List<String> sampleIds = entry.getValue();

            for (int i = 0; i < sampleIds.size(); ++i)
            {
                String sample1 = sampleIds.get(i);

                List<MultiBiopsyData> mbDataList1 = mSampleSvData.get(sample1);

                if (mbDataList1 == null || mbDataList1.isEmpty())
                    continue;

                for (int j = i + 1; j < sampleIds.size(); ++j)
                {
                    String sample2 = sampleIds.get(j);

                    List<MultiBiopsyData> mbDataList2 = mSampleSvData.get(sample2);

                    if (mbDataList2 == null || mbDataList2.isEmpty())
                        continue;

                    findSvMatches(mbDataList1, mbDataList2);
                }

                // cull any partials or non-exact if exact are found
                mbDataList1.forEach(x -> x.cullNonExactMatches());

                findSharedPrivateMerges(mbDataList1);

                for (int j = 0; j < sampleIds.size(); ++j)
                {
                    if(i == j)
                        continue;

                    String sample2 = sampleIds.get(j);
                    reportClusterOverlaps(sample1, mbDataList1, sample2);
                }
            }

            writeSampleSvData(patientId, sampleIds);
        }
    }

    private void findSvMatches(List<MultiBiopsyData> mbDataList1, List<MultiBiopsyData> mbDataList2)
    {
        // if this takes too long, use a reducible list
        for (MultiBiopsyData mbData1 : mbDataList1)
        {
            int partialMatches = 0;
            boolean hasSharedMatch = false;

            for (MultiBiopsyData mbData2 : mbDataList2)
            {
                String matchType = getMatchType(mbData1, mbData2);

                if (matchType == MATCH_TYPE_SHARED)
                {
                    mbData1.addMatchType(matchType, mbData2);
                    mbData2.addMatchType(matchType, mbData1);
                    hasSharedMatch = true;

                    // continue searching for shared matches (even though not expected)

                    // break;
                }
                else if (!hasSharedMatch && matchType == MATCH_TYPE_PARTIAL)
                {
                    mbData1.addMatchType(matchType, mbData2);
                    mbData2.addMatchType(matchType, mbData1);
                    ++partialMatches;

                    // continue searching within this sample
                }
            }

            if(partialMatches > 1 || (partialMatches > 0 && !mbData1.getSharedMatches().isEmpty()))
            {
                LOGGER.warn("sample({}) SV({}) found {} partial matches, {} shared",
                        mbData1.SampleId, mbData1.SvId, partialMatches, mbData1.getSharedMatches().size());
            }
        }
    }

    private String getMatchType(MultiBiopsyData mbData1, MultiBiopsyData mbData2)
    {
        // first check for an exact match on position and orientation
        boolean startMatched = areMatched(mbData1.ChrStart, mbData1.PosStart, mbData1.OrientStart,
                mbData2.ChrStart, mbData2.PosStart, mbData2.OrientStart);

        if (!startMatched)
        {
            // allow an SGL to match either end
            if (mbData1.Type == SGL && mbData2.Type != SGL)
            {
                if (areMatched(mbData1.ChrStart, mbData1.PosStart, mbData1.OrientStart, mbData2.ChrEnd, mbData2.PosEnd, mbData2.OrientEnd))
                    return MATCH_TYPE_PARTIAL;
            }
            else if (mbData2.Type == SGL && mbData1.Type != SGL)
            {
                if (areMatched(mbData2.ChrStart, mbData2.PosStart, mbData2.OrientStart, mbData1.ChrEnd, mbData1.PosEnd, mbData1.OrientEnd))
                    return MATCH_TYPE_PARTIAL;
            }

            return MATCH_TYPE_PRIVATE;
        }

        boolean endMatched = areMatched(mbData1.ChrEnd, mbData1.PosEnd, mbData1.OrientEnd,
                mbData2.ChrEnd, mbData2.PosEnd, mbData2.OrientEnd);

        if (endMatched)
            return MATCH_TYPE_SHARED;
        else
            return MATCH_TYPE_PRIVATE;
    }

    private void findSharedPrivateMerges(List<MultiBiopsyData> mbDataList)
    {
        // report details about any SV-pair in a cluster involving a shared with a private SV
        for(int i = 0; i < mbDataList.size(); ++i)
        {
            MultiBiopsyData mbData1 = mbDataList.get(i);

            for(int j = i+1; j < mbDataList.size(); ++j)
            {
                MultiBiopsyData mbData2 = mbDataList.get(j);

                if(mbData1.ClusterId != mbData2.ClusterId)
                    continue;

                if(mbData1.getMatchType() == mbData2.getMatchType())
                    continue;

                String reason = mbData1.getClusterReasonForId(mbData2.SvId);

                if(reason.isEmpty())
                    continue;

                writeSharedPrivateMerges(mbData1.SampleId, mbData1.ClusterId, mbData1.SvId, mbData1.getMatchType(),
                        mbData2.SvId, mbData2.getMatchType(), reason);
            }
        }
    }

    private void reportClusterOverlaps(String sample, List<MultiBiopsyData> mbDataList, String otherSample)
    {
        // for each cluster, report on the number of overlapping clusters from other same-patient samples
        // classify each cluster as: Private - all SVs are only in one sample, Exact - all shared SVs match,
        // SimpleSuperset - one sample has all the shared of a single other cluster, plus some private SVs
        // Subset - one sample has no private SVs, and some but not all the shared of a single other cluster
        // ComplexSuperset = one sample has no or some private, shared overlapping with more than 1 other cluster
        // otherwise Mixed

        int currentClusterId = -1;
        int currentClusterCount = 0;
        String currentResolvedType = "";
        List<Integer> matchingClusters = Lists.newArrayList();
        String otherClusterIds = "";
        int privateCount = 0;
        int sharedCount = 0;
        int otherClustersTotal = 0;

        for(int i = 0; i <= mbDataList.size(); ++i)
        {
            final MultiBiopsyData mbData = i < mbDataList.size() ? mbDataList.get(i) : null;

            if(mbData == null || currentClusterId != mbData.ClusterId)
            {
                if(i > 0)
                {
                    String overlapType;
                    if (sharedCount == 0)
                    {
                        overlapType = "Private";
                    }
                    else if (privateCount == 0 && matchingClusters.size() == 1 && currentClusterCount == otherClustersTotal)
                    {
                        overlapType = "Exact";
                    }
                    else if (privateCount == 0 && matchingClusters.size() == 1 && currentClusterCount < otherClustersTotal)
                    {
                        overlapType = "Subset";
                    }
                    else if (privateCount > 0 && matchingClusters.size() == 1 && sharedCount == otherClustersTotal)
                    {
                        overlapType = "SimpleSuperset";
                    }
                    else if (privateCount >= 0 && matchingClusters.size() > 1 && sharedCount == otherClustersTotal)
                    {
                        overlapType = "ComplexSuperset";
                    }
                    else
                    {
                        overlapType = "Mixed";
                    }

                    // write details to file
                    writeClusterOverlapData(sample, currentClusterId, currentClusterCount, currentResolvedType, privateCount, sharedCount,
                            matchingClusters.size(), otherClustersTotal, otherClusterIds, overlapType);
                }

                if(mbData == null)
                    break;

                matchingClusters.clear();
                currentClusterId = mbData.ClusterId;
                currentClusterCount= mbData.ClusterCount;
                currentResolvedType = mbData.ResolvedType;
                otherClusterIds = "";
                privateCount = 0;
                sharedCount = 0;
                otherClustersTotal = 0;
            }

            if(mbData.getMatchType() == MATCH_TYPE_PRIVATE)
            {
                ++privateCount;
                continue;
            }

            ++sharedCount;

            List<MultiBiopsyData> otherSvs = Lists.newArrayList(mbData.getSharedMatches());
            otherSvs.addAll(mbData.getPartialMatches());

            for (MultiBiopsyData otherSv : otherSvs)
            {
                if(!otherSv.SampleId.equals(otherSample))
                    continue;

                if(!matchingClusters.contains(otherSv.ClusterId))
                {
                    matchingClusters.add(otherSv.ClusterId);
                    otherClustersTotal += otherSv.ClusterCount;
                    otherClusterIds = appendStr(otherClusterIds, String.valueOf(otherSv.ClusterId), ';');
                }
            }
        }
    }

    private void writeSampleSvData(String patientId, List<String> sampleIds)
    {
        try
        {
            if (mSvWriter == null)
            {
                String outputFileName = mOutputDir + "LNX_MB_SV_DATA.csv";

                mSvWriter = createBufferedWriter(outputFileName, false);

                mSvWriter.write("PatientId,SampleId,SvId,Type,ClusterId,ClusterCount,ResolvedType,MatchType,OtherSvId");
                mSvWriter.newLine();
            }

            for(String sampleId : sampleIds)
            {
                List<MultiBiopsyData> mbDataList = mSampleSvData.get(sampleId);

                if(mbDataList == null | mbDataList.isEmpty())
                    continue;

                for(final MultiBiopsyData mbData : mbDataList)
                {
                    String outputStr = String.format("%s,%s,%d,%s,%d,%d,%s",
                            patientId, sampleId, mbData.SvId, mbData.Type,
                            mbData.ClusterId, mbData.ClusterCount, mbData.ResolvedType);

                    if(mbData.getMatchType() == MATCH_TYPE_PRIVATE)
                    {
                        mSvWriter.write(String.format("%s,%s,%d", outputStr, MATCH_TYPE_PRIVATE, -1));
                        mSvWriter.newLine();
                    }
                    else
                    {
                        for (MultiBiopsyData otherSv : mbData.getSharedMatches())
                        {
                            mSvWriter.write(String.format("%s,%s,%d", outputStr, MATCH_TYPE_SHARED, otherSv.SvId));
                            mSvWriter.newLine();
                        }

                        for (MultiBiopsyData otherSv : mbData.getPartialMatches())
                        {
                            mSvWriter.write(String.format("%s,%s,%d", outputStr, MATCH_TYPE_PARTIAL, otherSv.SvId));
                            mSvWriter.newLine();
                        }
                    }
                }
            }
        }
        catch (IOException e)
        {
            LOGGER.error("failed writing sv output data: {}", e.toString());
            return;
        }
    }

    private void writeSharedPrivateMerges(String sampleId, int clusterId, int svId1, String matchType1,
            int svId2, String matchType2, String clusterReason)
    {
        try
        {
            if (mMergeWriter == null)
            {
                String outputFileName = mOutputDir + "LNX_MB_MERGE_DATA.csv";

                mMergeWriter = createBufferedWriter(outputFileName, false);

                mMergeWriter.write("SampleId,ClusterId,SvId1,MatchType1,SvId2,MatchType2,ClusterReason");
                mMergeWriter.newLine();
            }

            mMergeWriter.write(String.format("%s,%d,%d,%s,%d,%s,%s",
                    sampleId, clusterId, svId1, matchType1, svId2, matchType2, clusterReason));
            mMergeWriter.newLine();
        }
        catch (IOException e)
        {
            LOGGER.error("failed writing merge output data: {}", e.toString());
            return;
        }
    }

    private void writeClusterOverlapData(final String sampleId, int clusterId, int clusterCount, final String resolvedType, int privateCount,
            int sharedCount, int matchingClustersCount, int otherClustersTotal, final String otherClusterIds, final String overlapType)
    {
        try
        {
            if (mClusterOverlapWriter == null)
            {
                String outputFileName = mOutputDir + "LNX_MB_CLUSTER_DATA.csv";

                mClusterOverlapWriter = createBufferedWriter(outputFileName, false);

                mClusterOverlapWriter.write("SampleId,ClusterId,ClusterCount,ResolvedType,PrivateCount,SharedCount");
                mClusterOverlapWriter.write(",MatchingClustersCount,OtherClustersTotal,OverlapType,OtherClusterIds");
                mClusterOverlapWriter.newLine();
            }

            mClusterOverlapWriter.write(String.format("%s,%d,%d,%s,%d,%d,%d,%d,%s,%s",
                    sampleId, clusterId, clusterCount, resolvedType, privateCount, sharedCount,
                    matchingClustersCount, otherClustersTotal, overlapType, otherClusterIds));

            mClusterOverlapWriter.newLine();
        }
        catch (IOException e)
        {
            LOGGER.error("failed writing merge output data: {}", e.toString());
            return;
        }
    }

    private void closeOutputFiles()
    {
        closeBufferedWriter(mSvWriter);
        closeBufferedWriter(mMergeWriter);
        closeBufferedWriter(mClusterOverlapWriter);
    }

}
