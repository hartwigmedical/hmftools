package com.hartwig.hmftools.statcalcs.cooc;

import static com.hartwig.hmftools.statcalcs.cooc.SampleCategoryData.SAMPLE_CAT_1_INDEX;
import static com.hartwig.hmftools.statcalcs.cooc.SampleCategoryData.SAMPLE_CAT_2_INDEX;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.stats.FisherExactTest;
import com.hartwig.hmftools.common.utils.FileWriterUtils;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// 3 variable generic co-occurence using Fisher's exact test
// first variable is a grouping variable (eg CancerType) and the other 2 variables are tested for co-occurence
// expected input: SampleId,Grouping,Category1,Category2
// eg SampleId,Gene,CancerType,LohType - would group by Gene and then test CancerType vs LohType correlation

public class ThreeVarCoOccurence
{
    private List<String> mSamples;

    // generic data structure for 3-way co-occurrence
    private final Map<String, List<SampleCategoryData>> mGroupingSampleGenericData;
    private final List<String> mGroupingValues;
    private final List<String> mCat1Values;
    private final List<String> mCat2Values;

    private String mGroupingField;
    private String mCategory1;
    private String mCategory2;

    private final FisherExactTest mFisherET;

    private BufferedWriter mWriter;

    private static final String THREE_VAR_INPUT_FILE = "three_var_input_file";

    private static final Logger LOGGER = LogManager.getLogger(ThreeVarCoOccurence.class);

    public ThreeVarCoOccurence(final CommandLine cmd, final String outputDir)
    {
        mGroupingField = "";
        mCategory1 = "";
        mCategory2 = "";

        final String inputFile = cmd.getOptionValue(THREE_VAR_INPUT_FILE);

        loadSampleGenericData(inputFile);

        final String outputFile = outputDir + "STATS_3VAR.csv";
        initialiseOutput(outputFile);

        mSamples = Lists.newArrayList();
        mGroupingSampleGenericData = new HashMap<>();
        mGroupingValues = Lists.newArrayList();
        mCat1Values = Lists.newArrayList();
        mCat2Values = Lists.newArrayList();

        mFisherET = new FisherExactTest();
    }

    public static void addCmdLineOptions(Options options)
    {
        options.addOption(THREE_VAR_INPUT_FILE, true, "Sample data with grouping and 2 variable");
    }

    public static boolean hasConfig(final CommandLine cmd)
    {
        return cmd.hasOption(THREE_VAR_INPUT_FILE);
    }

    private static String SPEC_GROUP_VAL = "";
    // private static String SPEC_GROUP_VAL = "SETD2";

    public void run()
    {
        if(mGroupingSampleGenericData.isEmpty())
            return;

        // for each of the group fields, calculate co-occurrence for each of the 2 categories
        mFisherET.initialise(mSamples.size());

        int hypothesesCount = mGroupingValues.size() * mCat1Values.size() * mCat2Values.size();

        for(final String groupingValue : mGroupingValues)
        {
            final List<SampleCategoryData> sampleDataList = mGroupingSampleGenericData.get(groupingValue);

            if(sampleDataList == null || sampleDataList.isEmpty())
                continue;

            int sampleCount = sampleDataList.size();

            LOGGER.info("processing group({}) with {} samples", groupingValue, sampleCount);

            if(groupingValue.equals(SPEC_GROUP_VAL))
            {
                LOGGER.debug("spec group value: {}", groupingValue);
            }

            for(final String cat1 : mCat1Values)
            {
                for(final String cat2 : mCat2Values)
                {
                    int withCat1 = 0;
                    int withCat2 = 0;
                    int withCat1WithCat2 = 0;
                    int withCat1NoCat2 = 0;
                    int noCat1WithCat2 = 0;
                    int noCat1NoCat2 = 0;

                    for(final SampleCategoryData sampleData : sampleDataList)
                    {
                        final List<String[]> catDataList = sampleData.getCategoryData();

                        boolean hasCat1 = false;
                        boolean hasCat2 = false;
                        boolean hasBoth = false;

                        for (final String[] catList : catDataList)
                        {
                            if (catList[SAMPLE_CAT_1_INDEX].equals(cat1) && catList[SAMPLE_CAT_2_INDEX].equals(cat2))
                            {
                                hasBoth = true;
                                break;
                            }

                            if (catList[SAMPLE_CAT_1_INDEX].equals(cat1))
                                hasCat1 = true;

                            if(catList[SAMPLE_CAT_2_INDEX].equals(cat2))
                                hasCat2 = true;
                        }

                        if(hasBoth)
                            ++withCat1WithCat2;
                        else if(hasCat1)
                            ++withCat1NoCat2;
                        else if(hasCat2)
                            ++noCat1WithCat2;
                        else
                            ++noCat1NoCat2;

                        if(hasCat1 || hasBoth)
                            ++withCat1;

                        if(hasCat2 || hasBoth)
                            ++withCat2;
                    }

                    double expectedVal = withCat1 / (double)sampleCount * withCat2;

                    double fisherProb = mFisherET.calc(withCat1WithCat2, noCat1WithCat2, withCat1NoCat2, noCat1NoCat2, expectedVal);

                    writeResultsData(groupingValue, cat1, cat2, sampleCount, withCat1, withCat2, fisherProb,
                            expectedVal, hypothesesCount, withCat1WithCat2, noCat1WithCat2, withCat1NoCat2, noCat1NoCat2);
                }
            }
        }

        FileWriterUtils.closeBufferedWriter(mWriter);
    }

    private boolean initialiseOutput(final String outputFileName)
    {
        try
        {
            mWriter = FileWriterUtils.createBufferedWriter(outputFileName, false);

            mWriter.write(String.format("%s,%s,%s,With%sCount", mGroupingField, mCategory1, mCategory2, mGroupingField));

            mWriter.write(String.format(",With%sCount,With%sCount,ExpectedCount,FETProb,TestCount,CountGtExp",
                    mCategory1, mCategory2));

            mWriter.write(String.format(",With%sWith%s,No%sWith%s,With%sNo%s,No%sNo%s",
                    mCategory1, mCategory2, mCategory1, mCategory2, mCategory1, mCategory2, mCategory1, mCategory2));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            LOGGER.error("failed to open stats output file: {}", e.toString());
            return false;
        }

        return true;
    }

    private void writeResultsData(final String groupingValue, final String cat1, final String cat2, int sampleCount,
            int withCat1, int withCat2, double fetProbability, double expectedVal, int testCount,
            int withCat1WithCat2, int noCat1WithCat2, int withCat1NoCat2, int noCat1NoCat2)
    {
        if (mWriter == null)
            return;

        try
        {
            mWriter.write(
                    String.format("%s,%s,%s,%d",
                            groupingValue, cat1, cat2, sampleCount));

            mWriter.write(
                    String.format(",%d,%d,%.2f,%4.3e,%d,%s,%d,%d,%d,%d",
                            withCat1, withCat2, expectedVal, fetProbability,
                            testCount, withCat1WithCat2 > expectedVal,
                            withCat1WithCat2, noCat1WithCat2, withCat1NoCat2, noCat1NoCat2));

            mWriter.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to stats output file: {}", e.toString());
        }
    }

    private static int GENERIC_DATA_CSV_COUNT = 4;
    private static int SAMPLE_INDEX = 0;
    private static int GROUPING_INDEX = 1;
    private static int CAT_1_INDEX = 2;
    private static int CAT_2_INDEX = 3;

    private void loadSampleGenericData(final String filename)
    {
        if (filename.isEmpty())
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            // field names are used for categories
            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("Empty driver genes CSV file({})", filename);
                return;
            }

            String[] items = line.split(",");

            if(items.length != GENERIC_DATA_CSV_COUNT)
            {
                LOGGER.error("invalid data line: {}", line);
                return;
            }

            mGroupingField = items[GROUPING_INDEX];
            mCategory1 = items[CAT_1_INDEX];
            mCategory2 = items[CAT_2_INDEX];

            int recordCount = 0;
            String currentGroupingValue = "";
            List<SampleCategoryData> sampleDataList = null;

            while ((line = fileReader.readLine()) != null)
            {
                // parse CSV data
                items = line.split(",");

                if (items.length != GENERIC_DATA_CSV_COUNT)
                    continue;

                ++recordCount;

                final String sampleId = items[SAMPLE_INDEX];
                final String groupingValue = items[GROUPING_INDEX];
                final String cat1Value = items[CAT_1_INDEX];
                final String cat2Value = items[CAT_2_INDEX];

                if(!currentGroupingValue.equals(groupingValue))
                {
                    currentGroupingValue = groupingValue;

                    sampleDataList = mGroupingSampleGenericData.get(groupingValue);

                    if(sampleDataList == null)
                    {
                        sampleDataList = Lists.newArrayList();
                        mGroupingSampleGenericData.put(groupingValue, sampleDataList);
                    }
                }

                if(!mSamples.contains(sampleId))
                    mSamples.add(sampleId);

                if(!mGroupingValues.contains(groupingValue))
                    mGroupingValues.add(groupingValue);

                if(!mCat1Values.contains(cat1Value))
                    mCat1Values.add(cat1Value);

                if(!mCat2Values.contains(cat2Value))
                    mCat2Values.add(cat2Value);

                boolean found = false;
                for(final SampleCategoryData sampleData : sampleDataList)
                {
                    if(sampleData.SampleId.equals(sampleId))
                    {
                        sampleData.addCategoryData(cat1Value, cat2Value);
                        found = true;
                        break;
                    }
                }

                if(!found)
                {
                    SampleCategoryData sampleData = new SampleCategoryData(sampleId);
                    sampleData.addCategoryData(cat1Value, cat2Value);
                    sampleDataList.add(sampleData);
                }
            }

            LOGGER.info("loaded {} sample generic data records", recordCount);

        }
        catch (IOException exception)
        {
            LOGGER.error("Failed to sample generic data file({})", filename);
        }
    }

}
