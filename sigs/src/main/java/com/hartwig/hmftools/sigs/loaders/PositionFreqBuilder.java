package com.hartwig.hmftools.sigs.loaders;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.sigs.PositionFrequencies.DEFAULT_POS_FREQ_MAX_SAMPLE_COUNT;
import static com.hartwig.hmftools.common.sigs.PositionFrequencies.getBucketIndex;
import static com.hartwig.hmftools.common.sigs.PositionFrequencies.getChromosomeFromIndex;
import static com.hartwig.hmftools.common.sigs.PositionFrequencies.getPositionFromIndex;
import static com.hartwig.hmftools.common.sigs.PositionFrequencies.initialisePositionCache;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.MatrixUtils.loadMatrixDataFile;
import static com.hartwig.hmftools.common.utils.MatrixUtils.writeMatrixData;
import static com.hartwig.hmftools.sigs.common.CommonUtils.LOG_DEBUG;
import static com.hartwig.hmftools.sigs.common.CommonUtils.OUTPUT_FILE_ID;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SIG_LOGGER;
import static com.hartwig.hmftools.sigs.common.CommonUtils.formOutputFilename;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class PositionFreqBuilder
{
    // position mapping
    private final Map<String,Integer> mChromosomeLengths;
    private final Map<String,Integer> mChromosomePosIndex;
    private final int mBucketSize;
    private final int mNewBucketSize;
    private int mPositionCacheSize;
    private final String mOutputDir;
    private final String mOutputFileId;
    private final String mSamplePosCountsFile;

    private final List<String> mSampleList;
    private final Map<String,int[]> mSamplePositionFrequencies;

    public static final String POSITION_BUCKET_SIZE = "position_bucket_size";
    public static final String CONVERT_POSITION_BUCKET_SIZE = "convert_pos_bucket_size";
    public static final String MAX_SAMPLE_COUNT = "max_sample_count";

    private static final String SAMPLE_DATA_FILE = "sample_data_file";
    private static final String POSITION_DATA_FILE = "position_data_file";

    public PositionFreqBuilder(final CommandLine cmd)
    {
        mChromosomeLengths = Maps.newHashMap();
        mChromosomePosIndex = Maps.newHashMap();
        mSamplePositionFrequencies = Maps.newHashMap();
        mSampleList = Lists.newArrayList();

        mBucketSize = Integer.parseInt(cmd.getOptionValue(POSITION_BUCKET_SIZE));
        mNewBucketSize = cmd.hasOption(CONVERT_POSITION_BUCKET_SIZE) ? Integer.parseInt(cmd.getOptionValue(CONVERT_POSITION_BUCKET_SIZE)) : 0;
        mPositionCacheSize = 0;
        mOutputDir = parseOutputDir(cmd);
        mOutputFileId = cmd.getOptionValue(OUTPUT_FILE_ID);

        mSamplePosCountsFile = cmd.getOptionValue(POSITION_DATA_FILE);

        loadSampleData(cmd.getOptionValue(SAMPLE_DATA_FILE));

        mPositionCacheSize = initialisePositionCache(mBucketSize, mChromosomeLengths, mChromosomePosIndex);

        SIG_LOGGER.info("position cache size({}) from position bucket position({})", mPositionCacheSize, mBucketSize);
    }

    public PositionFreqBuilder(final String outputDir, final String outputId, int positionBucketSize)
    {
        mChromosomeLengths = Maps.newHashMap();
        mChromosomePosIndex = Maps.newHashMap();
        mSamplePositionFrequencies = Maps.newHashMap();

        mBucketSize = positionBucketSize;
        mNewBucketSize = 0;
        mPositionCacheSize = 0;
        mOutputDir = outputDir;
        mOutputFileId = outputId;

        mPositionCacheSize = initialisePositionCache(mBucketSize, mChromosomeLengths, mChromosomePosIndex);
        SIG_LOGGER.info("position cache size({}) from position bucket position({})", mPositionCacheSize, mBucketSize);

        mSampleList = null;
        mSamplePosCountsFile = null;
    }


    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE_DATA_FILE, true, "File with sampleIds and cancer types");
        options.addOption(POSITION_BUCKET_SIZE, true, "Position bucket size");
        options.addOption(CONVERT_POSITION_BUCKET_SIZE, true, "Optional, convert from pos-freq matrix to larger bucket size");
        options.addOption(POSITION_DATA_FILE, true, "Position frequencies data file");
        options.addOption(MAX_SAMPLE_COUNT, true, "Max sample SNV count, default = 20K");
        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(OUTPUT_FILE_ID, true, "Output file ID");
        options.addOption(LOG_DEBUG, false, "Log verbose");
    }

    public void run()
    {
        if(mBucketSize <= 0 || mPositionCacheSize <= 0 || mSamplePosCountsFile == null)
        {
            SIG_LOGGER.error("failed to initialise");
            return;
        }

        if(mNewBucketSize > 0)
        {
            if(mNewBucketSize < mBucketSize)
            {
                SIG_LOGGER.warn("cannot shrink position frequencies from bucket size {} -> {}",
                        mBucketSize, mNewBucketSize);
                System.exit(1);
            }

            SIG_LOGGER.info("converting position frequencies from bucket size {} -> {}",
                    mBucketSize, mNewBucketSize);

            convertPosFrequencyBucketSize();
            return;
        }

        SIG_LOGGER.info("loading sample position counts");

        loadPositionFrequencies(mSamplePosCountsFile);

        if(mSamplePositionFrequencies.isEmpty())
        {
            SIG_LOGGER.error("no sample counts loaded");
            return;
        }

        SIG_LOGGER.info("writing sample position counts matrix");

        writeSampleData();

        SIG_LOGGER.info("position frequencies data sets complete");
    }

    private void convertPosFrequencyBucketSize()
    {
        final List<String> sampleNames = Lists.newArrayList();

        final Map<String,Integer> newChrPosIndexMap = Maps.newHashMap();
        int newPosCacheSize = initialisePositionCache(mNewBucketSize, mChromosomeLengths, newChrPosIndexMap);

        final Matrix oldSampleCounts = loadMatrixDataFile(mSamplePosCountsFile, sampleNames);
        final double[][] oldCounts = oldSampleCounts.getData();

        final Matrix newSampleCounts = new Matrix(newPosCacheSize, oldSampleCounts.Cols);
        final double[][] newCounts = newSampleCounts.getData();

        for(int r = 0; r < oldSampleCounts.Rows; ++r)
        {
            final String chromosome = getChromosomeFromIndex(mChromosomePosIndex, r);
            int position = getPositionFromIndex(mChromosomePosIndex, chromosome, r, mBucketSize);

            int newIndex = getBucketIndex(mNewBucketSize, newChrPosIndexMap, chromosome, position);

            for(int c = 0; c < oldSampleCounts.Cols; ++c)
            {
                newCounts[newIndex][c] += oldCounts[r][c];
            }
        }

        try
        {
            final String filename = String.format("%s/converted_pos_freq_counts_%d.csv", mOutputDir, mNewBucketSize);

            SIG_LOGGER.info("writing converted counts(samples={} posBuckets={}) to {}",
                    newSampleCounts.Cols, newSampleCounts.Rows, filename);

            BufferedWriter writer = createBufferedWriter(filename, false);

            writeMatrixData(writer, sampleNames, newSampleCounts, true);
            writer.close();
        }
        catch(IOException e)
        {
            SIG_LOGGER.error("error writing converted pos-freq counts: {}", e.toString());
        }
    }

    public void writeSampleCounts(final String sampleId, final Map<String,Map<Integer,Integer>> chrPositionCounts, boolean writeCoords)
    {
        SIG_LOGGER.info("sample({}) writing position frequencies data", sampleId);

        try
        {
            final String filename = mOutputDir + sampleId + ".sig.pos_freq_counts.csv";
            BufferedWriter writer = createBufferedWriter(filename, false);

            if(writeCoords)
                writer.write(String.format("Chromosome,Position,Frequency"));
            else
                writer.write(sampleId);

            writer.newLine();

            final int[] bucketFrequencies = new int[mPositionCacheSize];

            for(Map.Entry<String,Map<Integer,Integer>> chrEntry : chrPositionCounts.entrySet())
            {
                final String chromosome = chrEntry.getKey();

                for(Map.Entry<Integer, Integer> entry : chrEntry.getValue().entrySet())
                {
                    int position = entry.getKey();
                    int frequency = entry.getValue();

                    int posBucketIndex = getBucketIndex(mBucketSize, mChromosomePosIndex, chromosome, position);
                    bucketFrequencies[posBucketIndex] = frequency;
                }
            }

            for(int b = 0; b < mPositionCacheSize; ++b)
            {
                if(writeCoords)
                {
                    final String chromosome = getChromosomeFromIndex(mChromosomePosIndex, b);
                    int position = getPositionFromIndex(mChromosomePosIndex, chromosome, b, mBucketSize);

                    writer.write(String.format("%s,%d,%d", chromosome, position, bucketFrequencies[b]));
                }
                else
                {
                    writer.write(String.format("%d", bucketFrequencies[b]));
                }

                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            SIG_LOGGER.error("failed to write sample pos data output: {}", e.toString());
        }
    }

    private void loadSampleData(final String filename)
    {
        if(filename == null || filename.isEmpty())
            return;

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            if(fileData.isEmpty())
            {
                SIG_LOGGER.error("empty sample data file({})", filename);
                return;
            }

            final String header = fileData.get(0);
            fileData.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");

            for(final String line : fileData)
            {
                final String[] items = line.split(",", -1);
                String sampleId = items[fieldsIndexMap.get("SampleId")];

                mSampleList.add(sampleId);
            }

            SIG_LOGGER.info("loaded {} samples", mSampleList.size());
        }
        catch (IOException e)
        {
            SIG_LOGGER.error("failed to load sample data file({}): {}", filename, e.toString());
        }
    }

    private void loadPositionFrequencies(final String filename)
    {
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            // SampleId,Chromosome,Position,Count
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(line, ",");
            int sampleIndex = fieldsIndexMap.get("SampleId");
            int chrIndex = fieldsIndexMap.get("Chromosome");
            int posIndex = fieldsIndexMap.get("Position");
            int countIndex = fieldsIndexMap.get("Count");

            String currentSample = null;
            int[] currentCounts = null;

            line = fileReader.readLine(); // skip header

            int nextLog = 100;

            while (line != null)
            {
                // SampleId,SigName,SigContrib,SigPercent
                final String[] items = line.split(",", -1);
                String sampleId = items[sampleIndex];

                if(!mSampleList.contains(sampleId))
                {
                    line = fileReader.readLine();
                    continue;
                }

                String chromosome = items[chrIndex];
                int position = Integer.parseInt(items[posIndex]);
                int count = Integer.parseInt(items[countIndex]);

                if(currentSample == null || !currentSample.equals(sampleId))
                {
                    currentCounts = new int[mPositionCacheSize];
                    currentSample = sampleId;
                    mSamplePositionFrequencies.put(sampleId, currentCounts);

                    if(mSamplePositionFrequencies.size() >= nextLog)
                    {
                        SIG_LOGGER.debug("loaded {} sample position freq data", mSamplePositionFrequencies.size());
                        nextLog += 100;
                    }
                }

                int posBucketIndex = getBucketIndex(mBucketSize, mChromosomePosIndex, chromosome, position);
                currentCounts[posBucketIndex] += count;

                line = fileReader.readLine();
            }
        }
        catch (IOException e)
        {
            SIG_LOGGER.error("failed to read position frequencies file({}): {}", filename, e.toString());
        }
    }

    private void writeSampleData()
    {
        try
        {
            final String filename = formOutputFilename(mOutputDir, mOutputFileId, "pos_freq_matrix");

            BufferedWriter writer = createBufferedWriter(filename, false);

            final List<String> sampleIds = mSamplePositionFrequencies.keySet().stream().collect(Collectors.toList());
            writer.write(sampleIds.get(0));
            for(int i = 1; i < sampleIds.size(); ++i)
            {
                writer.write(String.format(",%s", sampleIds.get(i)));
            }

            writer.newLine();

            for(int b = 0; b < mPositionCacheSize; ++b)
            {
                writer.write(String.format("%d", mSamplePositionFrequencies.get(sampleIds.get(0))[b]));

                for(int i = 1; i < sampleIds.size(); ++i)
                {
                    writer.write(String.format(",%d", mSamplePositionFrequencies.get(sampleIds.get(i))[b]));
                }

                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            SIG_LOGGER.error("failed to write sample pos data output: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = new Options();
        PositionFreqBuilder.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        PositionFreqBuilder refDataBuilder = new PositionFreqBuilder(cmd);
        refDataBuilder.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
