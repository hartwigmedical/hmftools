package com.hartwig.hmftools.cup.utils;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.CuppaConfig.FLD_CANCER_TYPE;
import static com.hartwig.hmftools.cup.CuppaConfig.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.cup.CuppaConfig.REF_SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.cup.CuppaConfig.SUBSET_DELIM;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_BREAST;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_BREAST_TRIPLE_NEGATIVE;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_OTHER;
import static com.hartwig.hmftools.cup.utils.CompareUtils.EMPTY_RESULTS_CSV;
import static com.hartwig.hmftools.cup.utils.CompareUtils.resultInfoCsv;
import static com.hartwig.hmftools.cup.utils.CompareUtils.topRefResult;
import static com.hartwig.hmftools.cup.utils.CompareUtils.resultsMatch;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.cup.common.SampleResult;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class CuppaCompare
{
    private final Map<String,String> mOrigSampleCancerTypes;
    private final Map<String,String> mNewSampleCancerTypes;
    private final List<String> mClassifiers;
    private final String mFileOrig;
    private final String mFileNew;
    private final boolean mAllowCancerTypeRemap;

    private final Map<String,ResultCounts> mDataTypeResults;

    private final BufferedWriter mWriter;
    
    private static final String FILE_ORIG = "orig_file";
    private static final String FILE_NEW = "new_file";
    private static final String ORIG_SAMPLE_DATA = "orig_sample_data";
    private static final String NEW_SAMPLE_DATA = "new_sample_data";
    private static final String CLASSIFERS = "classifiers";
    private static final String ALLOW_CANCER_TYPE_REMAP = "allow_cancer_type_remap";

    public CuppaCompare(final CommandLine cmd)
    {
        mClassifiers = Lists.newArrayList();

        if(cmd.hasOption(CLASSIFERS))
        {
            String[] classifers = cmd.getOptionValue(CLASSIFERS).split(SUBSET_DELIM);
            Arrays.stream(classifers).forEach(x -> mClassifiers.add(x));
        }

        mFileOrig = cmd.getOptionValue(FILE_ORIG);
        mFileNew = cmd.getOptionValue(FILE_NEW);

        String outputDir = parseOutputDir(cmd);

        if(cmd.hasOption(REF_SAMPLE_DATA_FILE))
        {
            mOrigSampleCancerTypes = loadSampleData(cmd.getOptionValue(REF_SAMPLE_DATA_FILE));
            mNewSampleCancerTypes = mOrigSampleCancerTypes;
            mAllowCancerTypeRemap = false;
        }
        else
        {
            mOrigSampleCancerTypes = loadSampleData(cmd.getOptionValue(ORIG_SAMPLE_DATA));
            mNewSampleCancerTypes = loadSampleData(cmd.getOptionValue(NEW_SAMPLE_DATA));
            mAllowCancerTypeRemap = cmd.hasOption(ALLOW_CANCER_TYPE_REMAP);
        }

        mDataTypeResults = Maps.newHashMap();

        mWriter = initialiseWriter(outputDir, cmd.getOptionValue(OUTPUT_ID), mAllowCancerTypeRemap);
    }

    public void run()
    {
        if(mOrigSampleCancerTypes.isEmpty())
        {
            CUP_LOGGER.error("no samples loaded");
            System.exit(1);
        }

        if(mWriter == null)
        {
            System.exit(1);
            return;
        }

        if(mFileOrig == null || !Files.exists(Paths.get(mFileOrig)))
        {
            CUP_LOGGER.error("invalid original Cuppa file({})", mFileOrig);
            System.exit(1);
        }

        if(mFileNew == null || !Files.exists(Paths.get(mFileNew)))
        {
            CUP_LOGGER.error("invalid new Cuppa file({})", mFileNew);
            System.exit(1);
        }

        Map<String,List<SampleResult>> origAllResults = SampleResult.loadResults(mFileOrig);
        Map<String,List<SampleResult>> newAllResults = SampleResult.loadResults(mFileNew);

        if(origAllResults.isEmpty() || newAllResults.isEmpty())
        {
            CUP_LOGGER.warn("empty Cuppa results files: orig({}) new({})", origAllResults.size(), newAllResults.size());
            return;
        }

        Set<String> allSamples = mOrigSampleCancerTypes.keySet().stream().collect(Collectors.toSet());
        mNewSampleCancerTypes.keySet().forEach(x -> allSamples.add(x));

        CUP_LOGGER.info("loaded Cuppa data: orig(sample={} results={}) new(samples={} results={}) totalSamples({})",
                mOrigSampleCancerTypes.size(), origAllResults.values().stream().mapToInt(x -> x.size()).sum(),
                mNewSampleCancerTypes.size(), newAllResults.values().stream().mapToInt(x -> x.size()).sum(), allSamples.size());

        int processed = 0;

        for(String sampleId : allSamples)
        {
            if(mOrigSampleCancerTypes.containsKey(sampleId) && mNewSampleCancerTypes.containsKey(sampleId))
            {
                String origRefCancerType = mOrigSampleCancerTypes.get(sampleId);

                String newRefCancerType = mNewSampleCancerTypes.containsKey(sampleId) ? mNewSampleCancerTypes.get(sampleId) : origRefCancerType;

                if(!origRefCancerType.equals(newRefCancerType) && !mAllowCancerTypeRemap)
                {
                    CUP_LOGGER.warn("sample({}) has diff ref cancerTypes(orig={} new={}), skipping",
                            sampleId, origRefCancerType, newRefCancerType);
                    continue;
                }

                if(origRefCancerType.equals(CANCER_TYPE_OTHER))
                    continue;

                List<SampleResult> origResults = origAllResults.get(sampleId);
                List<SampleResult> newResults = newAllResults.get(sampleId);

                if(origResults == null || newResults == null)
                {
                    CUP_LOGGER.warn("sample({}) {} missing results, skipping", sampleId, origResults == null ? "original" : "new");
                    continue;
                }

                for(SampleResult origResult : origResults)
                {
                    SampleResult newResult = newResults.stream().filter(x -> resultsMatch(x, origResult)).findFirst().orElse(null);

                    processResults(sampleId, true, origRefCancerType, newRefCancerType, origResult, newResult);
                }
            }
            else
            {
                boolean hasOrig = mOrigSampleCancerTypes.containsKey(sampleId);

                List<SampleResult> results = hasOrig ? origAllResults.get(sampleId) : newAllResults.get(sampleId);

                if(results == null)
                {
                    CUP_LOGGER.warn("sample({}) {} missing results, skipping", sampleId, hasOrig ? "original" : "new");
                    continue;
                }

                String refCancerType = hasOrig ? mOrigSampleCancerTypes.get(sampleId) : mNewSampleCancerTypes.get(sampleId);

                if(refCancerType.equals(CANCER_TYPE_OTHER))
                    continue;

                for(SampleResult result : results)
                {
                    processResults(
                            sampleId, false, hasOrig ? refCancerType : "", !hasOrig ? refCancerType : "",
                            hasOrig ? result : null, !hasOrig ? result : null);
                }
            }

            ++processed;

            if((processed % 1000) == 0)
            {
                CUP_LOGGER.debug("compared {} samples", processed);
            }
        }

        for(Map.Entry<String,ResultCounts> entry : mDataTypeResults.entrySet())
        {
            String dataType = entry.getKey();
            ResultCounts resultCounts = entry.getValue();

            CUP_LOGGER.info(format("dataType(%s) bothCorrect(%d) bothIncorrect(%d) orig(correct=%d rate=%.3f) newCorrect(%d rate=%.3f)",
                    dataType, resultCounts.BothCorrect, resultCounts.BothIncorrect,
                    resultCounts.OrigCorrect, resultCounts.origRate(), resultCounts.NewCorrect, resultCounts.newRate()));
        }

        closeBufferedWriter(mWriter);
    }

    public static BufferedWriter initialiseWriter(final String outputDir, final String outputId, boolean writeNewCancerType)
    {
        try
        {
            String outputFileName = outputDir + "CUP_COMPARE";

            if(outputId != null)
                outputFileName += "." + outputId;

            outputFileName += ".csv";

            CUP_LOGGER.info("writing comparison results to {}", outputFileName);

            final BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("SampleId,Status,DataType,MatchType");

            if(writeNewCancerType)
                writer.write(",OrigRefCancerType,NewRefCancerType");
            else
                writer.write(",RefCancerType");

            writer.write(",OrigCancerType,OrigCancerValue,OrigRefCancerValue");
            writer.write(",NewCancerType,NewCancerValue,NewRefCancerValue");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to initialise Cuppa comparison file: {}", e.toString());
            return null;
        }
    }

    private static final String STATUS_BOTH = "BOTH";
    private static final String STATUS_ORIG_SAMPLE = "ORIG_SAMPLE_ONLY";
    private static final String STATUS_NEW_SAMPLE = "NEW_SAMPLE_ONLY";
    private static final String STATUS_ORIG_DATATYPE = "ORIG_DATATYPE_ONLY";
    private static final String STATUS_NEW_DATATYPE = "NEW_DATATYPE_ONLY";

    private static final String MATCH_TYPE_CORRECT = "CORRECT";
    private static final String MATCH_TYPE_ORIG_CORRECT = "ORIG_ONLY";
    private static final String MATCH_TYPE_NEW_CORRECT = "NEW_ONLY";
    private static final String MATCH_TYPE_INCORRECT = "INCORRECT";

    private boolean skipClassifier(final String dataType)
    {
        return !mClassifiers.isEmpty() && !mClassifiers.contains(dataType);
    }

    private void addResultCount(final String dataType, final String matchType)
    {
        ResultCounts results = mDataTypeResults.get(dataType);

        if(results == null)
        {
            results = new ResultCounts();
            mDataTypeResults.put(dataType, results);
        }

        if(matchType.equals(MATCH_TYPE_CORRECT))
            ++results.BothCorrect;
        else if(matchType.equals(MATCH_TYPE_INCORRECT))
            ++results.BothIncorrect;
        else if(matchType.equals(MATCH_TYPE_ORIG_CORRECT))
            ++results.OrigCorrect;
        else if(matchType.equals(MATCH_TYPE_NEW_CORRECT))
            ++results.NewCorrect;
    }

    private static boolean isCorrectSubtype(final String cancerType1, final String cancerType2)
    {
        return (cancerType1.equals(CANCER_TYPE_BREAST_TRIPLE_NEGATIVE) && cancerType2.equals(CANCER_TYPE_BREAST))
                || (cancerType1.equals(CANCER_TYPE_BREAST) && cancerType2.equals(CANCER_TYPE_BREAST_TRIPLE_NEGATIVE));
    }

    private void processResults(
            final String sampleId, boolean sampleInBoth, final String origRefCancerType, final String newRefCancerType,
            final SampleResult origResult, final SampleResult newResult)
    {
        try
        {
            if(origResult == null || newResult == null)
            {
                final SampleResult result = origResult != null ? origResult : newResult;
                final String refCancerType = origResult != null ? origRefCancerType : newRefCancerType;

                if(skipClassifier(result.DataType))
                    return;

                final String status = sampleInBoth ?
                        (origResult != null ? STATUS_ORIG_DATATYPE : STATUS_NEW_DATATYPE) :
                        (origResult != null ? STATUS_ORIG_SAMPLE : STATUS_NEW_SAMPLE);

                boolean isCorrect = topRefResult(result).equals(refCancerType) || isCorrectSubtype(topRefResult(result), refCancerType);
                String matchType = isCorrect ? MATCH_TYPE_CORRECT : MATCH_TYPE_INCORRECT;

                addResultCount(result.DataType, matchType);

                mWriter.write(format("%s,%s,%s,%s", sampleId, status, result.DataType, matchType));

                if(mAllowCancerTypeRemap)
                    mWriter.write(format(",%s,%s", origRefCancerType, newRefCancerType));
                else
                    mWriter.write(format(",%s", origRefCancerType));

                mWriter.write(format(",%s,%s",
                        origResult != null ? resultInfoCsv(origResult, refCancerType) : EMPTY_RESULTS_CSV,
                        newResult != null ? resultInfoCsv(newResult, refCancerType) : EMPTY_RESULTS_CSV));

                mWriter.newLine();
            }
            else
            {
                if(skipClassifier(origResult.DataType))
                    return;

                String topOrigType = topRefResult(origResult);
                String topNewType = topRefResult(newResult);
                boolean origCorrect = topOrigType.equals(origRefCancerType) || isCorrectSubtype(topOrigType, origRefCancerType);
                boolean newCorrect = topNewType.equals(newRefCancerType) || isCorrectSubtype(topNewType, newRefCancerType);

                String matchType;

                if(origCorrect && newCorrect)
                    matchType = MATCH_TYPE_CORRECT;
                else if(origCorrect)
                    matchType = MATCH_TYPE_ORIG_CORRECT;
                else if(newCorrect)
                    matchType = MATCH_TYPE_NEW_CORRECT;
                else
                    matchType = MATCH_TYPE_INCORRECT;

                addResultCount(origResult.DataType, matchType);

                mWriter.write(format("%s,%s,%s,%s", sampleId, STATUS_BOTH, origResult.DataType, matchType));

                if(mAllowCancerTypeRemap)
                    mWriter.write(format(",%s,%s", origRefCancerType, newRefCancerType));
                else
                    mWriter.write(format(",%s", origRefCancerType));

                mWriter.write(format(",%s,%s", resultInfoCsv(origResult, origRefCancerType), resultInfoCsv(newResult, newRefCancerType)));

                mWriter.newLine();
            }
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write Cuppa comparison file: {}", e.toString());
        }
    }

    private Map<String,String> loadSampleData(final String sampleDataFile)
    {
        final Map<String,String> sampleCancerTypes = Maps.newHashMap();

        try
        {
            final List<String> fileData = Files.readAllLines(new File(sampleDataFile).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);
            int sampleIdIndex = fieldsIndexMap.get(FLD_SAMPLE_ID);
            int cancerTypeIndex = fieldsIndexMap.get(FLD_CANCER_TYPE);

            for(final String line : fileData)
            {
                String[] values = line.split(DATA_DELIM, -1);
                sampleCancerTypes.put(values[sampleIdIndex], values[cancerTypeIndex]);
            }

            CUP_LOGGER.info("loaded {} samples from file: {}", sampleCancerTypes.size(), sampleDataFile);
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sample data file({}): {}", sampleDataFile, e.toString());
        }

        return sampleCancerTypes;
    }

    private class ResultCounts
    {
        public int BothCorrect;
        public int BothIncorrect;
        public int NewCorrect;
        public int OrigCorrect;

        public ResultCounts()
        {
            BothCorrect = 0;
            BothIncorrect = 0;
            NewCorrect = 0;
            OrigCorrect = 0;
        }

        public int total() { return BothCorrect + BothIncorrect + OrigCorrect + NewCorrect; }

        public double newRate() { return (BothCorrect + NewCorrect) / (double)total(); }
        public double origRate() { return (BothCorrect + OrigCorrect) / (double)total(); }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(CLASSIFERS, true, "Classifiers to compare");
        options.addOption(FILE_ORIG, true, "Original Cuppa results file");
        options.addOption(FILE_NEW, true, "New Cuppa results file");
        options.addOption(REF_SAMPLE_DATA_FILE, true, "Sample ref data file, if using the same for original and new results");
        options.addOption(ORIG_SAMPLE_DATA, true, "Original sample ref data file");
        options.addOption(NEW_SAMPLE_DATA, true, "New sample ref data file - optional if the same a original");
        options.addOption(ALLOW_CANCER_TYPE_REMAP, false, "Allow cancer types to change for separate sample ref data files");
        addLoggingOptions(options);
        addOutputOptions(options);

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        setLogLevel(cmd);

        CuppaCompare cuppaCompare = new CuppaCompare(cmd);
        cuppaCompare.run();
    }

}
