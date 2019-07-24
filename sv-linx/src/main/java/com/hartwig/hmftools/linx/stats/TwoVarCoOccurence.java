package com.hartwig.hmftools.linx.stats;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.stats.SampleCategoryData.SAMPLE_CAT_1_INDEX;
import static com.hartwig.hmftools.linx.stats.SampleCategoryData.SAMPLE_CAT_2_INDEX;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// 2 variable generic co-occurence using Fisher's exact test
// the 2 variables are tested for co-occurence
// expected input: SampleId,Category1,Category2,Count
// eg SampleId,Gene,LohType - would test CancerType vs LohType correlation

public class TwoVarCoOccurence
{
    private List<String> mSamples;

    // generic data structure for 2-way co-occurrence
    private List<TwoCategoryData> mCategoryCountsData;
    private List<String> mCat1Values;
    private List<String> mCat2Values;

    private String mCategory1;
    private String mCategory2;

    private FisherExactTest mFisherET;

    private BufferedWriter mWriter;

    private static final Logger LOGGER = LogManager.getLogger(ThreeVarCoOccurence.class);

    public TwoVarCoOccurence()
    {
        mCategory1 = "";
        mCategory2 = "";

        mSamples = Lists.newArrayList();
        mCategoryCountsData = Lists.newArrayList();
        mCat1Values = Lists.newArrayList();
        mCat2Values = Lists.newArrayList();

        mFisherET = new FisherExactTest();

        mWriter = null;
    }

    public void run()
    {
        if(mCategoryCountsData.isEmpty())
            return;

        int totalRecords = mCategoryCountsData.stream().mapToInt(x -> x.Count).sum();

        mFisherET.initialise(totalRecords);

        int hypothesesCount = mCat1Values.size() * mCat2Values.size();

        LOGGER.info("processing {} 2-var records", totalRecords);

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

        closeBufferedWriter(mWriter);
    }

    private boolean initialiseOutput(final String outputFileName)
    {
        try
        {
            mWriter = createBufferedWriter(outputFileName, false);

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

    public boolean initialise(final String filename, final String outputDir)
    {
        loadSampleGenericData(filename);

        final String outputFile = outputDir + "SVA_STATS_2VAR.csv";
        return initialiseOutput(outputFile);
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
                LOGGER.error("empy 2-var CSV file({})", filename);
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
