package com.hartwig.hmftools.isofox.novel.cohort;

import static java.lang.Math.log;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.ALT_SPLICE_JUNCTION;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunction.FLD_ALT_SJ_FRAG_COUNT;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunction.FLD_ALT_SJ_POS_END;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunction.FLD_ALT_SJ_POS_START;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunction.FLD_ALT_SJ_TYPE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_CHROMOSOME;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_ID;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class AltSjCohortMatrix
{
    private final CohortConfig mConfig;

    // other config
    private final int mMinSampleThreshold;
    private final int mMinFragments;

    private final Map<String,Map<Integer,Integer>> mAltSpliceJunctionSampleData; // map of alt-SJ to map of sampleIndex to frag counts
    private final Map<String,AltSjLocation> mAltSpliceJunctionData;
    private final Map<Integer,String> mMatrixIndexMap;
    private final Map<String,Integer> mMatrixKeyMap;
    private final Map<String,Integer> mItemMatrixIndexMap; // map of sample or cancer-type into matrix
    private final Map<String,Integer> mSampleIndexMap;

    private final boolean mGroupByCancer;
    private final List<String> mCancerTypes;
    private final int mExpectMatrixSize;
    private int[][] mMatrixData;

    private static final String ALT_SJ_BY_CANCER = "alt_sj_by_cancer";
    private static final String ALT_SJ_MIN_SAMPLES = "alt_sj_min_samples";
    private static final String ALT_SJ_MIN_FRAGS = "alt_sj_min_frags";
    private static final String ALT_SJ_EXP_COUNT = "alt_sj_expected_count";

    public AltSjCohortMatrix(final CohortConfig config, final CommandLine cmd)
    {
        mConfig = config;
        mAltSpliceJunctionSampleData = Maps.newHashMap();
        mAltSpliceJunctionData = Maps.newHashMap();
        mMatrixIndexMap = Maps.newHashMap();
        mMatrixKeyMap = Maps.newHashMap();
        mSampleIndexMap = Maps.newHashMap();

        mExpectMatrixSize = Integer.parseInt(cmd.getOptionValue(ALT_SJ_EXP_COUNT, "10000"));

        mGroupByCancer = cmd.hasOption(ALT_SJ_BY_CANCER);

        mItemMatrixIndexMap = Maps.newHashMap();
        mCancerTypes = mConfig.SampleData.CancerTypeSamples.keySet().stream().collect(Collectors.toList());

        if(mGroupByCancer)
        {
            for(int s = 0; s < mConfig.SampleData.SampleIds.size(); ++s)
            {
                final String sampleId = mConfig.SampleData.SampleIds.get(s);
                final String cancerType = mConfig.SampleData.SampleCancerType.get(sampleId);

                for(int c = 0; c < mCancerTypes.size(); ++c)
                {
                    if(cancerType.equals(mCancerTypes.get(c)))
                    {
                        mItemMatrixIndexMap.put(sampleId, c);
                        break;
                    }
                }
            }
        }

        for(int s = 0; s < mConfig.SampleData.SampleIds.size(); ++s)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(s);
            mSampleIndexMap.put(sampleId, s);

            if(!mGroupByCancer)
                mItemMatrixIndexMap.put(sampleId, s);
        }

        int itemCount = mGroupByCancer ? mCancerTypes.size() : mConfig.SampleData.SampleIds.size();
        mMatrixData = new int[mExpectMatrixSize][itemCount];

        mMinSampleThreshold = Integer.parseInt(cmd.getOptionValue(ALT_SJ_MIN_SAMPLES, "1"));
        mMinFragments = Integer.parseInt(cmd.getOptionValue(ALT_SJ_MIN_FRAGS, "1"));
    }

    public static void addCmdLineOptions(final Options options)
    {
        options.addOption(ALT_SJ_MIN_SAMPLES, true, "Min number of samples to report an alt SJ");
        options.addOption(ALT_SJ_MIN_FRAGS, true, "Min frag count supporting alt-SJs");
        options.addOption(ALT_SJ_EXP_COUNT, true, "Expected count of common alt-SJs");
        options.addOption(ALT_SJ_BY_CANCER, false, "Group by cancer type");
    }

    private static final int SHIFT_CHECK = 10;
    private static final int LOG_CHECK = 100;

    public void processAltSpliceJunctions()
    {
        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, ALT_SPLICE_JUNCTION, filenames))
            return;

        int nextLog = LOG_CHECK;
        int nextShift = SHIFT_CHECK;

        final Map<String,Integer> fieldsMap = Maps.newHashMap();

        // load each sample's alt SJs and consolidate into a single list
        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path altSJFile = filenames.get(i);

            loadFile(sampleId, altSJFile, fieldsMap);

            if(i >= nextShift)
            {
                nextShift += SHIFT_CHECK;
                shiftMapDataToMatrix(1000);
            }

            if(i >= nextLog)
            {
                nextLog += LOG_CHECK;
                ISF_LOGGER.info("processed {} samples, matrixCount({}) mapCount({})",
                        i, mMatrixIndexMap.size(), mAltSpliceJunctionSampleData.size());
            }
        }

        shiftMapDataToMatrix(0);

        ISF_LOGGER.info("writing matrix data for {} altSJs and {} {}",
                mMatrixIndexMap.size(), mGroupByCancer ? mCancerTypes.size() : mConfig.SampleData.SampleIds.size(),
                mGroupByCancer ? "cancer types" : "samples");

        writeMatrixData();
    }

    public void loadFile(final String sampleId, final Path filename, final Map<String,Integer> fieldsIndexMap)
    {
        try
        {
            final List<String> lines = Files.readAllLines(filename);

            if(fieldsIndexMap.isEmpty())
                fieldsIndexMap.putAll(createFieldsIndexMap(lines.get(0), DELIMITER));

            int geneIndex = fieldsIndexMap.get(FLD_GENE_ID);
            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posStartIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_START);
            int posEndIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_END);
            int typeIndex = fieldsIndexMap.get(FLD_ALT_SJ_TYPE);
            int fragCountIndex = fieldsIndexMap.get(FLD_ALT_SJ_FRAG_COUNT);

            Integer itemIndex = mItemMatrixIndexMap.get(sampleId);
            Integer sampleIndex = mSampleIndexMap.get(sampleId);

            int matrixAdded = 0;
            int skipped = 0;
            int mapAdded = 0;

            for(int i = 1; i < lines.size(); ++i)
            {
                final String[] items = lines.get(i).split(DELIMITER);
                final AltSjLocation altSJ = AltSjLocation.fromCsv(items, geneIndex, chrIndex, posStartIndex, posEndIndex, typeIndex);
                final int fragCount = Integer.parseInt(items[fragCountIndex]);

                Integer matrixIndex = mMatrixKeyMap.get(altSJ.Key);

                if(matrixIndex != null)
                {
                    mMatrixData[matrixIndex][itemIndex] += fragCount;
                    ++matrixAdded;
                }
                else
                {
                    Map<Integer,Integer> altSampleData = mAltSpliceJunctionSampleData.get(altSJ.Key);

                    if(altSampleData == null)
                    {
                        // skip the first alt-SJ if less than min fragments
                        if(fragCount < mMinFragments)
                        {
                            ++skipped;
                            continue;
                        }

                        altSampleData = Maps.newHashMap();
                        mAltSpliceJunctionSampleData.put(altSJ.Key, altSampleData);
                        mAltSpliceJunctionData.put(altSJ.Key, altSJ);
                    }

                    altSampleData.put(sampleIndex, fragCount);
                    ++mapAdded;
                }
            }

            ISF_LOGGER.debug("sample({}) loaded {} alt-SJ records, skipped({}) matrix({}) map({})",
                    sampleId, lines.size() - 1, skipped, matrixAdded, mapAdded);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to alt splice junction load file({}): {}", filename.toString(), e.toString());
            return;
        }
    }

    private void shiftMapDataToMatrix(int minShiftCount)
    {
        final List<String> shiftKeys = mAltSpliceJunctionSampleData.entrySet().stream()
                .filter(x -> x.getValue().size() >= mMinSampleThreshold)
                .map(x -> x.getKey())
                .collect(Collectors.toList());

        if(shiftKeys.size() < minShiftCount)
            return;

        int shiftCount = shiftKeys.size();
        int itemCount = mGroupByCancer ? mCancerTypes.size() : mConfig.SampleData.SampleIds.size();
        int exitingMatrixRowCount = mMatrixIndexMap.size();
        int newMatrixRowCount = exitingMatrixRowCount + shiftCount;

        if(newMatrixRowCount > mMatrixData.length)
        {
            // grow the matrix and copy in the new data
            int newExpectedSize = max(newMatrixRowCount, mMatrixData.length + mExpectMatrixSize);
            ISF_LOGGER.debug("growing matrix({} -> {})", mMatrixData.length, newExpectedSize);

            final int[][] newMatrix = new int[newExpectedSize][itemCount];

            for(int i = 0; i < exitingMatrixRowCount; ++i)
            {
                for(int s = 0; s < itemCount; ++s)
                {
                    newMatrix[i][s] = mMatrixData[i][s];
                }
            }

            mMatrixData = newMatrix;
        }

        for(String asjKey : shiftKeys)
        {
            int asjIndex = mMatrixIndexMap.size();
            Map<Integer,Integer> sampleData = mAltSpliceJunctionSampleData.get(asjKey);

            for(Map.Entry<Integer,Integer> entry : sampleData.entrySet())
            {
                Integer sampleIndex = entry.getKey();
                String sampleId = mConfig.SampleData.SampleIds.get(sampleIndex);
                Integer itemIndex = mItemMatrixIndexMap.get(sampleId);
                mMatrixData[asjIndex][itemIndex] += entry.getValue();
            }

            mMatrixIndexMap.put(asjIndex, asjKey);
            mMatrixKeyMap.put(asjKey, asjIndex);
            mAltSpliceJunctionSampleData.remove(asjKey);
        }

        ISF_LOGGER.debug("shifted {} asjSJs to matrix, matrixSize({}) mapSize({})",
                shiftKeys.size(), mMatrixIndexMap.size(),  mAltSpliceJunctionData.size());
    }

    private void writeMatrixData()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("alt_sj_cohort_matrix.csv");
            final BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("GeneId,Chromosome,PosStart,PosEnd");

            final List<String> itemNames = mGroupByCancer ? mCancerTypes : mConfig.SampleData.SampleIds;

            for(String sampleId : itemNames)
            {
                writer.write(String.format(",%s", sampleId));
            }

            writer.newLine();

            // write matrix data to file
            for(int asjIndex = 0; asjIndex < mMatrixIndexMap.size(); ++asjIndex)
            {
                final String asjKey = mMatrixIndexMap.get(asjIndex);
                final AltSjLocation altSJ = mAltSpliceJunctionData.get(asjKey);

                writer.write(String.format("%s,%s,%d,%d",
                        altSJ.GeneId, altSJ.Location.Chromosome, altSJ.Location.start(), altSJ.Location.end()));

                for(int s = 0; s < itemNames.size(); ++s)
                {
                    int fragCount = mMatrixData[asjIndex][s];
                    writer.write(String.format(",%d", fragCount));
                }

                writer.newLine();
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write alt-SJ cohort matrix file: {}", e.toString());
        }
    }
}
