package com.hartwig.hmftools.isofox.novel.cohort;

import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_FRAG_COUNT;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_END;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_START;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.ALT_SPLICE_JUNCTION;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
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
    // private final int mMinSampleThreshold;
    private final int mMinFragments;

    private final Map<String,AltSjLocation> mAltSjDataMap;
    private final List<AltSjLocation> mAltSjDataList;
    private final Map<String,Integer> mAltSjMatrixIndexMap; // map of alt-SJ key to matrix index (ie row)

    private final Map<String,Integer> mSampleMatrixIndexMap; // map of sample into matrix
    private final Map<String,Integer> mCancerMatrixIndexMap; // map of cancer-type into matrix

    private final List<String> mCancerTypes;
    private int[][] mSampleMatrixData;
    private int[][] mCancerMatrixData;

    private static final String ALT_SJ_COHORT_SITES_FILE = "alt_sj_cohort_sites_file";
    private static final String ALT_SJ_WRITE_CANCER_MATRIX = "alt_sj_write_cancer_matrix";
    private static final String ALT_SJ_WRITE_SAMPLE_MATRIX = "alt_sj_write_sample_matrix";
    private static final String ALT_SJ_MIN_FRAGS = "alt_sj_min_frags";

    public AltSjCohortMatrix(final CohortConfig config, final CommandLine cmd)
    {
        mConfig = config;
        mAltSjDataMap = Maps.newHashMap();
        mAltSjMatrixIndexMap = Maps.newHashMap();
        mAltSjDataList = Lists.newArrayList();

        mMinFragments = Integer.parseInt(cmd.getOptionValue(ALT_SJ_MIN_FRAGS, "1"));

        if(cmd.hasOption(ALT_SJ_COHORT_SITES_FILE))
        {
            loadCohortSites(cmd.getOptionValue(ALT_SJ_COHORT_SITES_FILE));
        }

        int altSjSiteCount = mAltSjDataMap.size();

        mCancerMatrixIndexMap = Maps.newHashMap();
        mCancerMatrixData = null;

        mCancerTypes = mConfig.SampleData.CancerTypeSamples.keySet().stream().collect(Collectors.toList());

        if(cmd.hasOption(ALT_SJ_WRITE_CANCER_MATRIX))
        {
            int cancerCount = mCancerTypes.size();
            mCancerMatrixData = new int[altSjSiteCount][cancerCount];

            for(int i = 0; i < mCancerTypes.size(); ++i)
            {
                String cancerType = mCancerTypes.get(i);
                mCancerMatrixIndexMap.put(cancerType, i);
            }
        }

        mSampleMatrixIndexMap = Maps.newHashMap();
        mSampleMatrixData = null;

        if(cmd.hasOption(ALT_SJ_WRITE_SAMPLE_MATRIX))
        {
            int sampleCount = mConfig.SampleData.SampleIds.size();
            mSampleMatrixData = new int[altSjSiteCount][sampleCount];

            for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
            {
                final String sampleId = mConfig.SampleData.SampleIds.get(i);
                mSampleMatrixIndexMap.put(sampleId, i);
            }
        }
    }

    public static void addCmdLineOptions(final Options options)
    {
        options.addOption(ALT_SJ_COHORT_SITES_FILE, true, "Alt-SJ reoccurring sites in cohort to filter by");
        options.addOption(ALT_SJ_WRITE_CANCER_MATRIX, false, "Write cancer matrix for cohort alt-SJs");
        options.addOption(ALT_SJ_WRITE_SAMPLE_MATRIX, false, "Write sample matrix for cohort alt-SJs");
        options.addOption(ALT_SJ_MIN_FRAGS, true, "Min frag count supporting alt-SJs");
    }

    private static final int LOG_CHECK = 100;

    public void processAltSpliceJunctions()
    {
        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, ALT_SPLICE_JUNCTION, filenames))
            return;

        int nextLog = LOG_CHECK;

        // load each sample's alt SJs and consolidate into a single list
        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path altSJFile = filenames.get(i);

            loadFile(sampleId, altSJFile);

            if(i >= nextLog)
            {
                nextLog += LOG_CHECK;
                ISF_LOGGER.info("processed {} samples", i);
            }
        }

        if(mCancerMatrixData != null)
        {
            ISF_LOGGER.info("writing matrix data for {} cancer types", mCancerTypes.size());
            writeMatrixData(mCancerMatrixData, mCancerTypes, "cancer");
        }

        if(mSampleMatrixData != null)
        {
            ISF_LOGGER.info("writing matrix data for {} samples", mConfig.SampleData.SampleIds.size());
            writeMatrixData(mSampleMatrixData, mConfig.SampleData.SampleIds, "sample");
        }
    }

    public void loadCohortSites(final String filename)
    {
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String header = fileReader.readLine();
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIMITER);

            int geneIndex = fieldsIndexMap.get(FLD_GENE_ID);
            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posStartIndex = fieldsIndexMap.get("PosStart");
            int posEndIndex = fieldsIndexMap.get("PosEnd");

            String line = fileReader.readLine();

            while(line != null)
            {
                final String[] items = line.split(DELIMITER, -1);

                final AltSjLocation altSJ = AltSjLocation.fromCsv(items, geneIndex, chrIndex, posStartIndex, posEndIndex);

                mAltSjDataMap.put(altSJ.Key, altSJ);
                mAltSjDataList.add(altSJ);
                mAltSjMatrixIndexMap.put(altSJ.Key, mAltSjMatrixIndexMap.size());
                line = fileReader.readLine();
            }

            ISF_LOGGER.info("loaded {} cohort alt-SJ sites from file", mAltSjDataList.size(), filename);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to cohort alt-SJ site file({}): {}", filename.toString(), e.toString());
            return;
        }
    }

    public void loadFile(final String sampleId, final Path filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(filename);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);

            int geneIndex = fieldsIndexMap.get(FLD_GENE_ID);
            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posStartIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_START);
            int posEndIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_END);
            int fragCountIndex = fieldsIndexMap.get(FLD_ALT_SJ_FRAG_COUNT);

            Integer sampleMatrixIndex = mSampleMatrixIndexMap.get(sampleId);
            String cancerType = mConfig.SampleData.SampleCancerType.get(sampleId);
            Integer cancerMatrixIndex = mCancerMatrixIndexMap.get(cancerType);

            int matrixAdded = 0;

            for(int i = 1; i < lines.size(); ++i)
            {
                final String[] items = lines.get(i).split(DELIMITER,-1);

                AltSjLocation altSJ = null;

                try
                {
                    altSJ = AltSjLocation.fromCsv(items, geneIndex, chrIndex, posStartIndex, posEndIndex);
                }
                catch(Exception e)
                {
                    ISF_LOGGER.error("sample({}) invalid alt-SJ line: {}", sampleId, lines.get(i));
                    return;
                }

                Integer asjMatrixIndex = mAltSjMatrixIndexMap.get(altSJ.Key);

                if(asjMatrixIndex == null)
                    continue;

                final int fragCount = Integer.parseInt(items[fragCountIndex]);

                if(fragCount < mMinFragments)
                    continue;

                ++matrixAdded;

                if(sampleMatrixIndex != null)
                {
                    mSampleMatrixData[asjMatrixIndex][sampleMatrixIndex] += fragCount;
                }

                if(cancerMatrixIndex != null)
                {
                    mCancerMatrixData[asjMatrixIndex][cancerMatrixIndex] += fragCount;
                }
            }

            ISF_LOGGER.debug("sample({}) loaded {} alt-SJ records, added({})",
                    sampleId, lines.size() - 1, matrixAdded);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to alt splice junction load file({}): {}", filename.toString(), e.toString());
            return;
        }
    }

    private void writeMatrixData(final int[][] matrixData, final List<String> itemNames, final String fileId)
    {
        try
        {
            final String fileName = String.format("alt_sj_cohort_matrix_%s.csv", fileId);
            final String outputFileName = mConfig.formCohortFilename(fileName);
            final BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("GeneId,Chromosome,PosStart,PosEnd");

            for(String sampleId : itemNames)
            {
                writer.write(String.format(",%s", sampleId));
            }

            writer.newLine();

            // write matrix data to file
            for(int asjIndex = 0; asjIndex < mAltSjDataList.size(); ++asjIndex)
            {
                final AltSjLocation altSJ = mAltSjDataList.get(asjIndex);

                writer.write(String.format("%s,%s,%d,%d",
                        altSJ.GeneId, altSJ.Location.Chromosome, altSJ.Location.start(), altSJ.Location.end()));

                for(int s = 0; s < itemNames.size(); ++s)
                {
                    int fragCount = matrixData[asjIndex][s];
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
