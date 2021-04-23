package com.hartwig.hmftools.statcalcs.cooc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.stats.FisherExactTest;
import com.hartwig.hmftools.common.utils.FileWriterUtils;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// 2 variable generic co-occurence using Fisher's exact test
// the 2 variables are tested for co-occurence
// expected input: Category1,Category2,Count
// eg Gene,LohType,Count - would test Gene vs LohType correlation

public class TwoVarCoOccurence
{
    // generic data structure for 2-way co-occurrence
    private final List<TwoCategoryData> mCategoryCountsData;
    private final List<String> mCat1Values;
    private final List<String> mCat2Values;

    private String mCategory1;
    private String mCategory2;

    private final FisherExactTest mFisherET;

    private BufferedWriter mWriter;

    private static final String TWO_VAR_INPUT_FILE = "two_var_input_file";

    private static final Logger LOGGER = LogManager.getLogger(TwoVarCoOccurence.class);

    public TwoVarCoOccurence(final CommandLine cmd, final String outputDir)
    {
        mCategory1 = "";
        mCategory2 = "";

        mCategoryCountsData = Lists.newArrayList();
        mCat1Values = Lists.newArrayList();
        mCat2Values = Lists.newArrayList();

        final String inputFile = cmd.getOptionValue(TWO_VAR_INPUT_FILE);
        loadSampleGenericData(inputFile);

        mFisherET = new FisherExactTest();

        final String outputFile = outputDir + "STATS_2VAR.csv";
        initialiseOutput(outputFile);
    }

    public static void addCmdLineOptions(Options options)
    {
        options.addOption(TWO_VAR_INPUT_FILE, true, "Sample data with 2 variables");
    }

    public static boolean hasConfig(final CommandLine cmd)
    {
        return cmd.hasOption(TWO_VAR_INPUT_FILE);
    }

    public void run()
    {
        if(mCategoryCountsData.isEmpty())
            return;

        int totalRecords = mCategoryCountsData.stream().mapToInt(x -> x.Count).sum();

        mFisherET.initialise(totalRecords);

        int hypothesesCount = mCat1Values.size() * mCat2Values.size();

        LOGGER.info("processing {} 2-var records, hypothese({})", totalRecords, hypothesesCount);

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

                for(final TwoCategoryData catData : mCategoryCountsData)
                {
                    boolean hasCat1 = false;
                    boolean hasCat2 = false;
                    boolean hasBoth = false;

                    if (catData.Category1.equals(cat1) && catData.Category2.equals(cat2))
                    {
                        withCat1WithCat2 += catData.Count;
                        hasBoth = true;
                    }
                    else if (catData.Category1.equals(cat1))
                    {
                        withCat1NoCat2 += catData.Count;
                        hasCat1 = true;
                    }
                    else if(catData.Category2.equals(cat2))
                    {
                        noCat1WithCat2 += catData.Count;
                        hasCat2 = true;
                    }

                    if(hasCat1 || hasBoth)
                        withCat1 += catData.Count;

                    if(hasCat2 || hasBoth)
                        withCat2 += catData.Count;
                }

                noCat1NoCat2 = totalRecords - withCat1WithCat2 - noCat1WithCat2 - withCat1NoCat2;

                double expectedVal = withCat1 / (double)totalRecords * withCat2;

                double fisherProb = mFisherET.calc(withCat1WithCat2, noCat1WithCat2, withCat1NoCat2, noCat1NoCat2, expectedVal);

                writeResultsData(cat1, cat2, totalRecords, withCat1, withCat2, fisherProb,
                        expectedVal, hypothesesCount, withCat1WithCat2, noCat1WithCat2, withCat1NoCat2, noCat1NoCat2);
            }
        }

        FileWriterUtils.closeBufferedWriter(mWriter);
    }

    private boolean initialiseOutput(final String outputFileName)
    {
        try
        {
            mWriter = FileWriterUtils.createBufferedWriter(outputFileName, false);

            mWriter.write(String.format("%s,%s,TotalCount", mCategory1, mCategory2));

            mWriter.write(String.format(",With%s,With%s,ExpectedCount,FETProb,TestCount,CountGtExp",
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

    private void writeResultsData(final String cat1, final String cat2, int totalCount,
            int withCat1, int withCat2, double fetProbability, double expectedVal, int testCount,
            int withCat1WithCat2, int noCat1WithCat2, int withCat1NoCat2, int noCat1NoCat2)
    {
        if (mWriter == null)
            return;

        try
        {
            mWriter.write(String.format("%s,%s,%d", cat1, cat2, totalCount));

            mWriter.write(String.format(",%d,%d,%.2f,%4.3e,%d,%s,%d,%d,%d,%d",
                    withCat1, withCat2, expectedVal, fetProbability, testCount, withCat1WithCat2 > expectedVal,
                    withCat1WithCat2, noCat1WithCat2, withCat1NoCat2, noCat1NoCat2));

            mWriter.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to stats output file: {}", e.toString());
        }
    }

    private static int CAT_1_INDEX = 0;
    private static int CAT_2_INDEX = 1;
    private static int COUNT_INDEX = 2;
    private static int GENERIC_DATA_CSV_COUNT = COUNT_INDEX + 1;

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
                LOGGER.error("empty 2-var CSV file({})", filename);
                return;
            }

            String[] items = line.split(",");

            if(items.length != GENERIC_DATA_CSV_COUNT)
            {
                LOGGER.error("invalid data line: {}", line);
                return;
            }

            mCategory1 = items[CAT_1_INDEX];
            mCategory2 = items[CAT_2_INDEX];

            int recordCount = 0;

            while ((line = fileReader.readLine()) != null)
            {
                // parse CSV data
                items = line.split(",");

                if (items.length != GENERIC_DATA_CSV_COUNT)
                    continue;

                ++recordCount;

                final String cat1Value = items[CAT_1_INDEX];
                final String cat2Value = items[CAT_2_INDEX];
                int counts = Integer.parseInt(items[COUNT_INDEX]);

                if(!mCat1Values.contains(cat1Value))
                    mCat1Values.add(cat1Value);

                if(!mCat2Values.contains(cat2Value))
                    mCat2Values.add(cat2Value);

                mCategoryCountsData.add(new TwoCategoryData(cat1Value, cat2Value, counts));
            }

            LOGGER.info("loaded {} 2-var data records", recordCount);

        }
        catch (IOException exception)
        {
            LOGGER.error("failed to load 2-var data file({})", filename);
        }
    }

}
