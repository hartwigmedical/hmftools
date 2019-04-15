package com.hartwig.hmftools.svanalysis.stats;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.SvAnalyser.DATA_OUTPUT_PATH;
import static com.hartwig.hmftools.svanalysis.types.MultiBiopsyData.MATCH_TYPE_PARTIAL;
import static com.hartwig.hmftools.svanalysis.types.MultiBiopsyData.MATCH_TYPE_PRIVATE;
import static com.hartwig.hmftools.svanalysis.types.MultiBiopsyData.MATCH_TYPE_SHARED;

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
import com.hartwig.hmftools.svanalysis.types.MultiBiopsyData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class MultipleBiopsyAnalyser
{
    private Map<String, List<String>> mPatientSampleIdsMap;
    private Map<String, List<MultiBiopsyData>> mSampleSvData;

    private BufferedWriter mSvWriter;
    private BufferedWriter mMergeWriter;
    private String mOutputDir;

    private static final Logger LOGGER = LogManager.getLogger(MultipleBiopsyAnalyser.class);

    public MultipleBiopsyAnalyser()
    {
        mPatientSampleIdsMap = Maps.newHashMap();
        mSampleSvData = Maps.newHashMap();
        mSvWriter = null;
        mMergeWriter = null;
        mOutputDir = "";
    }

    private static String PATIENT_SAMPLE_IDS_FILE = "patient_ids_file";
    private static String SVA_INPUT_FILE = "sva_svs_file";

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(PATIENT_SAMPLE_IDS_FILE, true, "File mapping PatientIds to SampleIds file");
        options.addOption(SVA_INPUT_FILE, true, "SVA SVs file");
    }

    public boolean loadData(final CommandLine cmd)
    {
        if (!loadPatientSampleData(cmd.getOptionValue(PATIENT_SAMPLE_IDS_FILE)))
            return false;

        if (!loadSampleSVData(cmd.getOptionValue(SVA_INPUT_FILE)))
            return false;

        mOutputDir = cmd.getOptionValue(DATA_OUTPUT_PATH);

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
        } catch (IOException e)
        {
            LOGGER.error("Failed to load patient sample IDs file({}): {}", filename, e.toString());
        }

        return true;
    }

    public void runAnalysis()
    {
        mapSamples();
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

                findSharedPrivateMerges(mbDataList1);
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

            for (MultiBiopsyData mbData2 : mbDataList2)
            {
                String matchType = getMatchType(mbData1, mbData2);

                if (matchType == MATCH_TYPE_SHARED)
                {
                    mbData1.addMatchType(matchType, mbData2.SvId);
                    mbData2.addMatchType(matchType, mbData1.SvId);
                    break;
                }
                else if (matchType == MATCH_TYPE_PARTIAL)
                {
                    mbData1.addMatchType(matchType, mbData2.SvId);
                    mbData2.addMatchType(matchType, mbData1.SvId);
                    ++partialMatches;

                    // continue searching within this sample
                }
            }

            if(partialMatches > 1 || (partialMatches > 0 && !mbData1.getSharedMatches().isEmpty()))
            {
                LOGGER.warn("Sample({}) SV({}) found {} partial matches, {} shared",
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

    private boolean areMatched(String chr1, long pos1, byte orient1, String chr2, long pos2, byte orient2)
    {
        return areMatched(chr1, pos1, orient1, chr2, pos2, orient2, false);
    }

    private static final int MAX_POS_DIFF = 5;

    private boolean areMatched(String chr1, long pos1, byte orient1, String chr2, long pos2, byte orient2, boolean requireExact)
    {
        if (!chr1.equals(chr2) || orient1 != orient2)
            return false;

        long posDiff = abs(pos1 - pos2);

        if (requireExact)
            return posDiff == 0;
        else
            return posDiff <= MAX_POS_DIFF;
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

    private void writeSampleSvData(String patientId, List<String> sampleIds)
    {
        try
        {
            if (mSvWriter == null)
            {
                String outputFileName = mOutputDir + "SVA_MB_SV_DATA.csv";

                mSvWriter = createBufferedWriter(outputFileName, false);

                mSvWriter.write("PatientId,SampleId,SvId,Type,ClusterCount,ResolvedType,MatchType,OtherSvId");
                mSvWriter.newLine();
            }

            for(String sampleId : sampleIds)
            {
                List<MultiBiopsyData> mbDataList = mSampleSvData.get(sampleId);

                if(mbDataList == null | mbDataList.isEmpty())
                    continue;

                for(final MultiBiopsyData mbData : mbDataList)
                {
                    String outputStr = String.format("%s,%s,%d,%s,%d,%s",
                            patientId, sampleId, mbData.SvId, mbData.Type,
                            mbData.ClusterCount, mbData.ResolvedType);

                    if(mbData.getMatchType() == MATCH_TYPE_PRIVATE)
                    {
                        mSvWriter.write(String.format("%s,%s,%d", outputStr, MATCH_TYPE_PRIVATE, -1));
                        mSvWriter.newLine();
                    }
                    else
                    {
                        for (Integer otherSvId : mbData.getSharedMatches())
                        {
                            mSvWriter.write(String.format("%s,%s,%d", outputStr, MATCH_TYPE_SHARED, otherSvId));
                            mSvWriter.newLine();
                        }

                        for (Integer otherSvId : mbData.getPartialMatches())
                        {
                            mSvWriter.write(String.format("%s,%s,%d", outputStr, MATCH_TYPE_PARTIAL, otherSvId));
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
                String outputFileName = mOutputDir + "SVA_MB_MERGE_DATA.csv";

                mMergeWriter = createBufferedWriter(outputFileName, false);

                mMergeWriter.write("SampleId,ClusterId,SvId1,MatchType1,SvId2Type,MatchType2,ClusterReason");
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

}
