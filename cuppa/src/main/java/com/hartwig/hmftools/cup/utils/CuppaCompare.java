package com.hartwig.hmftools.cup.utils;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.cup.CuppaConfig.CANCER_SUBTYPE_OTHER;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.CuppaConfig.FLD_CANCER_TYPE;
import static com.hartwig.hmftools.cup.CuppaConfig.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.cup.CuppaConfig.REF_SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.cup.CuppaConfig.SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_OTHER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.SampleResult;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class CuppaCompare
{
    private final Map<String,String> mSampleCancerTypes;
    private final String mFileOrig;
    private final String mFileNew;

    private final BufferedWriter mWriter;
    
    private static final String FILE_ORIG = "orig_file";
    private static final String FILE_NEW = "new_file";

    public CuppaCompare(final CommandLine cmd)
    {
        mSampleCancerTypes = Maps.newHashMap();
        mFileOrig = cmd.getOptionValue(FILE_ORIG);
        mFileNew = cmd.getOptionValue(FILE_NEW);

        String outputDir = parseOutputDir(cmd);

        mWriter = initialiseWriter(outputDir, cmd.getOptionValue(OUTPUT_ID));

        loadRefSampleData(cmd.getOptionValue(REF_SAMPLE_DATA_FILE));
    }

    public void run()
    {
        if(mSampleCancerTypes.isEmpty())
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

        CUP_LOGGER.info("loaded Cuppa results orig({}) new({})", origAllResults.size(), newAllResults.size());

        int matchedCount = 0;
        int diffCount = 0;
        int compared = 0;

        int unmatchedSamples = 0;
        int unmatchedResults = 0;
        int bothCorrect = 0;
        int origOnlyCorrect = 0;
        int newOnlyCorrect = 0;
        int bothIncorrect = 0;
        int missingCancerTypes = 0;

        for(Map.Entry<String,String> entry : mSampleCancerTypes.entrySet())
        {
            String sampleId = entry.getKey();
            String refCancerType = entry.getValue();

            if(refCancerType.equals(CANCER_TYPE_OTHER))
                continue;

            List<SampleResult> origResults = origAllResults.get(sampleId);
            List<SampleResult> newResults = newAllResults.get(sampleId);

            if(newResults == null) // skip samples not present in both cohort results
            {
                ++unmatchedSamples;
                continue;
            }

            for(SampleResult origResult : origResults)
            {
                SampleResult newResult = newResults.stream().filter(x -> x.typeMatches(origResult)).findFirst().orElse(null);

                if(newResult == null)
                {
                    ++unmatchedResults;
                    continue;
                }

                String origCancerType = origResult.topRefResult();
                String newCancerType = newResult.topRefResult();

                boolean origCorrect = origCancerType.equals(refCancerType);
                boolean newCorrect = newCancerType.equals(refCancerType);

                if(origCorrect && newCorrect)
                {
                    ++bothCorrect;
                    continue;
                }

                if(origCorrect)
                    ++origOnlyCorrect;
                else if(newCorrect)
                    ++newOnlyCorrect;
                else
                    ++bothIncorrect;

                if(!origResult.CancerTypeValues.containsKey(refCancerType) || !newResult.CancerTypeValues.containsKey(refCancerType))
                {
                    ++missingCancerTypes;
                    continue;
                }

                writeDiffs(
                        sampleId, refCancerType, origResult.Category, origResult.DataType,
                        origCancerType, origResult.CancerTypeValues.get(origCancerType), origResult.CancerTypeValues.get(refCancerType),
                        newCancerType, newResult.CancerTypeValues.get(newCancerType), newResult.CancerTypeValues.get(refCancerType));
            }

            ++compared;

            if((compared % 1000) == 0)
            {
                CUP_LOGGER.debug("compared {} samples", compared);
            }
        }

        CUP_LOGGER.info("bothCorrect({}) origCorrect({}) newCorrect({}) neither({}) unmatched(samples={} cancerTypes={})",
                bothCorrect, origOnlyCorrect, newOnlyCorrect, bothIncorrect, unmatchedSamples, missingCancerTypes);

        closeBufferedWriter(mWriter);
    }

    public static BufferedWriter initialiseWriter(final String outputDir, final String outputId)
    {
        try
        {
            String outputFileName = outputDir + "CUP_COMPARE";

            if(outputId != null)
                outputFileName += "." + outputId;

            outputFileName += ".csv";

            final BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("SampleId,Category,DataType,RefCancerType,IsCorrect");
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

    private void writeDiffs(
            final String sampleId, final String refCancerType, final CategoryType category, final String dataType,
            final String origCancerType, double origCancerTypeValue, double origRefCancerTypeValue,
            final String newCancerType, double newCancerTypeValue, double newRefCancerTypeValue)
    {
        try
        {
            boolean origCorrect = origCancerType.equals(refCancerType);
            boolean newCorrect = newCancerType.equals(refCancerType);
            String result = origCorrect && newCorrect ? "BOTH" : (origCorrect ? "ORIG" : (newCorrect ? "NEW" : "NEITHER"));
            mWriter.write(String.format("%s,%s,%s,%s,%s", sampleId, category, dataType, refCancerType, result));
            mWriter.write(String.format(",%s,%s,%s", origCancerType, origCancerTypeValue, origRefCancerTypeValue));
            mWriter.write(String.format(",%s,%s,%s", newCancerType, newCancerTypeValue, newRefCancerTypeValue));
            mWriter.newLine();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to initialise Cuppa comparison file: {}", e.toString());
        }
    }

    private void loadRefSampleData(final String sampleDataFile)
    {
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
                mSampleCancerTypes.put(values[sampleIdIndex], values[cancerTypeIndex]);
            }

            CUP_LOGGER.info("loaded {} sample and cancer types from file: {}", mSampleCancerTypes.size(), sampleDataFile);
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sample data file({}): {}", sampleDataFile, e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(FILE_ORIG, true, "Original Cuppa results file");
        options.addOption(FILE_NEW, true, "New Cuppa results file");
        options.addOption(REF_SAMPLE_DATA_FILE, true, "Sample ref data file");
        addLoggingOptions(options);
        addOutputOptions(options);

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        setLogLevel(cmd);

        CuppaCompare cuppaCompare = new CuppaCompare(cmd);
        cuppaCompare.run();
    }

}
