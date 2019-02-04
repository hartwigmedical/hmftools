package com.hartwig.hmftools.svanalysis.stats;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class StatisticRoutines
{
    private List<DriveGeneData> mDriveGeneData;
    private List<SampleCountsData> mSampleCountsData;

    private List<String> mSamples;
    private List<String> mGenes;
    private List<String> mCancerTypes;
    private List<String> mCategories;

    private Map<String, List<SampleData>> mCancerSampleData;
    private Map<String, Integer> mCategoryIndexMap;
    private Map<String, Integer> mGeneIndexMap;

    private int[][][] mSampleCountsMatrix;

    private FisherExactTest mFisherET;

    private BufferedWriter mWriter;

    private static final Logger LOGGER = LogManager.getLogger(StatisticRoutines.class);

    public StatisticRoutines()
    {
        mDriveGeneData = Lists.newArrayList();
        mSampleCountsData = Lists.newArrayList();
        mCancerTypes = Lists.newArrayList();
        mSamples = Lists.newArrayList();
        mGenes = Lists.newArrayList();
        mCategories = Lists.newArrayList();
        mCategoryIndexMap = new HashMap();
        mGeneIndexMap = new HashMap();
        mCancerSampleData = new HashMap();
        mSampleCountsMatrix = null;
        mFisherET = new FisherExactTest();
        mWriter = null;
    }

    private static String DRIVER_GENES_FILE = "driver_genes_file";
    private static String SAMPLE_COUNTS_FILE = "sample_counts_file";
    private static String OUTPUT_FILE = "stats_results_file";

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(DRIVER_GENES_FILE, true, "Drive genes file");
        options.addOption(SAMPLE_COUNTS_FILE, true, "Sample counts file");
        options.addOption(OUTPUT_FILE, true, "Results file");
    }

    public boolean loadConfig(final CommandLine cmd)
    {
        loadDriverGeneData(cmd.getOptionValue(DRIVER_GENES_FILE));
        loadSampleCountsData(cmd.getOptionValue(SAMPLE_COUNTS_FILE));

        return initialiseOutput(cmd.getOptionValue(OUTPUT_FILE));
    }

    private static String VALUE_TRUE = "TRUE";
    private static String VALUE_FALSE = "FALSE";
    private static String VALUE_UNCLEAR = "UNCLEAR";

    private static int ENR_CAT_NO_GENE = 0;
    private static int UNC_CAT_NO_GENE = 1;
    private static int ENR_CAT_WITH_GENE = 2;
    private static int UNC_CAT_WITH_GENE = 3;
    private static int ENR_CAT_UNC_GENE = 4;
    private static int UNC_CAT_UNC_GENE = 5;
    private static int INDEX_COUNT = 6;

    private static String SPEC_CANCER = "";
    // private static String SPEC_CANCER = "Eye";
    private static String SPEC_GENE = "";
    // private static String SPEC_GENE = "OR11H1";
    private static String SPEC_CATEGORY = "";
    // private static String SPEC_CATEGORY = "DUP_LT_100";

    public void runStatistics()
    {


        List<SampleData> allSampleDataList = Lists.newArrayList();

        for (final String cancerType : mCancerTypes)
        {
            if (!SPEC_CANCER.isEmpty() && !cancerType.equals(SPEC_CANCER))
                continue;

            List<SampleData> sampleDataList = mCancerSampleData.get(cancerType);
            allSampleDataList.addAll(sampleDataList);

            analyseCancerType(cancerType, sampleDataList);

        }
        
        // repeat for all cancers combined
        analyseCancerType("All", allSampleDataList);

        closeBufferedWriter(mWriter);

        LOGGER.info("analysis complete");
    }

    private void analyseCancerType(final String cancerType, final List<SampleData> sampleDataList)
    {
        /* calc method (for each cancer type):
        - each sample data has enriched categories and driver genes
        - from this it can contribute to each non-false cell
        - also tally up across all samples the counts per category and gene
         */
        int categoryCount = mCategories.size();
        int geneCount = mGenes.size();

        int sampleCount = sampleDataList.size();
        LOGGER.info("processing cancerType({}) with {} samples", cancerType, sampleCount);

        // populate the matrix with counts
        mSampleCountsMatrix = new int[geneCount][categoryCount][INDEX_COUNT];

        mFisherET.initialise(sampleCount);

        Map<String, Integer> withCategoryTotals = new HashMap();
        Map<String, Integer> unclearCategoryTotals = new HashMap();
        Map<String, Integer> withGeneTotals = new HashMap();
        Map<String, Integer> unclearGeneTotals = new HashMap();

        for(final SampleData sampleData : sampleDataList)
        {
            // calc category totals
            for(int i = 0; i <= 1; ++i)
            {
                final List<String> categories = (i == 0) ? sampleData.CategoryKnown : sampleData.CategoryUnclear;
                Map<String, Integer> totalsMap = (i == 0) ? withCategoryTotals : unclearCategoryTotals;

                for (final String category : categories)
                {
                    Integer total = totalsMap.get(category);

                    if(total == null)
                        totalsMap.put(category, 1);
                    else
                        totalsMap.put(category, total+1);
                }
            }

            // and then gene totals
            for(int i = 0; i <= 1; ++i)
            {
                final List<String> genes = (i == 0) ? sampleData.GeneKnown : sampleData.GeneUnclear;
                Map<String, Integer> totalsMap = (i == 0) ? withGeneTotals : unclearGeneTotals;

                for (final String gene : genes)
                {
                    Integer total = totalsMap.get(gene);

                    if(total == null)
                        totalsMap.put(gene, 1);
                    else
                        totalsMap.put(gene, total+1);
                }
            }

            for(int i = 0; i <= 1; ++i)
            {
                final List<String> categories = (i == 0) ? sampleData.CategoryKnown : sampleData.CategoryUnclear;

                for (int j = 0; j <= 1; ++j)
                {
                    final List<String> genes = (j == 0) ? sampleData.GeneKnown : sampleData.GeneUnclear;

                    int countsIndex;
                    if(i == 0 && j == 0)
                        countsIndex = ENR_CAT_WITH_GENE;
                    else if(i == 0 && j == 1)
                        countsIndex = ENR_CAT_UNC_GENE;
                    else if(i == 1 && j == 0)
                        countsIndex = UNC_CAT_WITH_GENE;
                    else
                        countsIndex = UNC_CAT_UNC_GENE;

                    for(final String category : categories)
                    {
                        for(final String gene : genes)
                        {
                            int categoryIndex = mCategoryIndexMap.get(category);
                            int geneIndex = mGeneIndexMap.get(gene);

                            ++mSampleCountsMatrix[geneIndex][categoryIndex][countsIndex];
                        }
                    }
                }
            }
        }

        LOGGER.info("cancerType({}) input counts populated", cancerType);

        for(final String category : mCategories)
        {
            if(!SPEC_CATEGORY.isEmpty() && !cancerType.equals(SPEC_CATEGORY))
                continue;

            int withCatTotal = withCategoryTotals.containsKey(category) ? withCategoryTotals.get(category) : 0;
            int uncCatTotal = unclearCategoryTotals.containsKey(category) ? unclearCategoryTotals.get(category) : 0;
            int noCatTotal = sampleCount - withCatTotal - uncCatTotal;

            int categoryIndex = mCategoryIndexMap.get(category);

            for(final String gene : mGenes)
            {
                if(!SPEC_GENE.isEmpty() && !cancerType.equals(SPEC_GENE))
                    continue;

                int withGeneTotal = withGeneTotals.containsKey(gene) ? withGeneTotals.get(gene) : 0;
                int uncGeneTotal = unclearGeneTotals.containsKey(gene) ? unclearGeneTotals.get(gene) : 0;
                int noGeneTotal = sampleCount - withGeneTotal - uncGeneTotal;

                int geneIndex = mGeneIndexMap.get(gene);

                int withCatWithGene = mSampleCountsMatrix[geneIndex][categoryIndex][ENR_CAT_WITH_GENE];
                int withCatUncGene = mSampleCountsMatrix[geneIndex][categoryIndex][ENR_CAT_UNC_GENE];
                int uncCatWithGene = mSampleCountsMatrix[geneIndex][categoryIndex][UNC_CAT_WITH_GENE];
                int uncCatUncGene = mSampleCountsMatrix[geneIndex][categoryIndex][UNC_CAT_UNC_GENE];

                // infer the others
                int noCatWithGene = withGeneTotal - withCatWithGene - uncCatWithGene;
                int noCatUncGene = uncGeneTotal - withCatUncGene - uncCatUncGene;
                int withCatNoGene = withCatTotal - withCatWithGene - withCatUncGene;
                int uncCatNoGene = uncCatTotal - uncCatWithGene - uncCatUncGene;
                int noCatNoGene = noGeneTotal - withCatNoGene - uncCatNoGene;

                if(withCatWithGene < 0 || noCatWithGene < 0 || withCatNoGene < 0 || noCatNoGene < 0)
                {
                    LOGGER.warn("INVALID COUNTS: cancer({}) samples({}) gene({}) counts(w={} u={} n={}) cat({}) counts(w={}) u={} n={})",
                            cancerType, sampleCount, gene, withGeneTotal, uncGeneTotal, noGeneTotal,
                            category, withCatTotal, uncCatTotal, noCatTotal);

                    LOGGER.warn("with cat: total({}) withGene({}) uncGene({}) noGene({})",
                            withCatTotal, withCatWithGene, withCatUncGene, withCatNoGene);

                    LOGGER.warn("unclear cat: total({}) withGene({}) uncGene({}) noGene({})",
                            uncCatTotal, uncCatWithGene, uncCatUncGene, uncCatNoGene);

                    LOGGER.warn("no cat: total({}) withGene({}) uncGene({}) noGene({})",
                            noCatTotal, noCatWithGene, noCatUncGene, noCatNoGene);

                    return;
                }

                /*
                if(gene.contains("BRCA"))
                {
                    LOGGER.debug("gene({}) counts(w={} u={} n={}) cat({}) counts(w={}) u={} n={})",
                            gene, withGeneTotal, uncGeneTotal, noGeneTotal,
                            category, withCatTotal, uncCatTotal, noCatTotal);

                    LOGGER.debug("with cat: total({}) withGene({}) uncGene({}) noGene({})",
                            withCatTotal, withCatWithGene, withCatUncGene, withCatNoGene);

                    LOGGER.debug("unclear cat: total({}) withGene({}) uncGene({}) noGene({})",
                            uncCatTotal, uncCatWithGene, uncCatUncGene, uncCatNoGene);

                    LOGGER.debug("no cat: total({}) withGene({}) uncGene({}) noGene({})",
                            noCatTotal, noCatWithGene, noCatUncGene, noCatNoGene);
                }
                */

                double geneSamplesPerc = withGeneTotal/(double)sampleCount;
                double expectedVal  = withCatTotal * geneSamplesPerc;
                double fisherProb = calcFisherExact(withCatWithGene, noCatWithGene, withCatNoGene, noCatNoGene, expectedVal);

                writeResultsData(cancerType, gene, category, sampleCount, withGeneTotal, withCatTotal, fisherProb,
                        expectedVal, withCatWithGene, noCatWithGene, withCatNoGene, noCatNoGene);
            }
        }

        LOGGER.info("cancerType({}) results written to file", cancerType);
    }

    private boolean initialiseOutput(final String outputFileName)
    {
        try
        {
            mWriter = createBufferedWriter(outputFileName, false);

            mWriter.write("CancerType,Gene,Category,SampleCount");
            mWriter.write(",WithGeneCount,WithCategoryCount,ExpectedCount,FETProb");
            mWriter.write(",WithCatWithGene,NoCatWithGene,WithCatNoGene,NoCatNoGene");

            mWriter.newLine();
        }
        catch(IOException e)
        {
            LOGGER.error("failed to open stats output file: {}", e.toString());
            return false;
        }

        return true;
    }

    private void writeResultsData(final String cancerType, final String gene, final String category, int sampleCount,
            int withGeneTotal, int withCatTotal, double fetProbability, double expectedVal,
            int withCatWithGene, int noCatWithGene, int withCatNoGene, int noCatNoGene)
    {
        if (mWriter == null)
            return;

        try
        {
            mWriter.write(
                    String.format("%s,%s,%s,%d",
                            cancerType, gene, category, sampleCount));

            mWriter.write(
                    String.format(",%d,%d,%.2f,%4.3e,%d,%d,%d,%d",
                            withGeneTotal, withCatTotal, expectedVal, fetProbability,
                            withCatWithGene, noCatWithGene, withCatNoGene, noCatNoGene));

            mWriter.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to stats output file: {}", e.toString());
        }
    }


    private double calcFisherExact(int withAwithB, int withANoB, int noAWithB, int noAnoB, double expectedCount)
    {
        if(withAwithB > expectedCount)
        {
            return mFisherET.getRightTailedP(withAwithB, noAWithB, withANoB, noAnoB);
        }
        else
        {
            return mFisherET.getLeftTailedP(withAwithB, noAWithB, withANoB, noAnoB);
        }
    }

    private SampleData getOrCreateSampleData(final String cancerType, final String sampleId)
    {
        List<SampleData> sampleDataList = mCancerSampleData.get(cancerType);

        SampleData sampleData = null;
        if (sampleDataList == null)
        {
            sampleDataList = Lists.newArrayList();
            mCancerSampleData.put(cancerType, sampleDataList);
        }
        else
        {
            for (final SampleData sd : sampleDataList)
            {
                if (sd.SampleId.equals(sampleId))
                {
                    sampleData = sd;
                    break;
                }
            }
        }

        if (sampleData == null)
        {
            sampleData = new SampleData(sampleId);
            sampleDataList.add(sampleData);
        }

        return sampleData;
    }

    private void loadDriverGeneData(final String filename)
    {
        if (filename.isEmpty())
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            // skip field names
            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("Empty driver genes CSV file({})", filename);
                return;
            }

            SampleData currentSampleData = null;

            while ((line = fileReader.readLine()) != null)
            {
                // parse CSV data
                String[] items = line.split(",");

                if (items.length != 4)
                    continue;

                DriveGeneData data = new DriveGeneData(
                        items[0], items[1], items[2], items[3]);

                mDriveGeneData.add(data);

                if(!mCancerTypes.contains(data.CancerType))
                    mCancerTypes.add(data.CancerType);

                if(!mSamples.contains(data.SampleId))
                    mSamples.add(data.SampleId);

                if(!mGenes.contains(data.Gene))
                {
                    mGeneIndexMap.put(data.Gene, mGenes.size());
                    mGenes.add(data.Gene);
                }

                if(currentSampleData == null || !currentSampleData.SampleId.equals(data.SampleId))
                {
                    currentSampleData = getOrCreateSampleData(data.CancerType, data.SampleId);
                }

                if(data.DriverStatus.equals(VALUE_TRUE))
                    currentSampleData.GeneKnown.add(data.Gene);
                else if(data.DriverStatus.equals(VALUE_UNCLEAR))
                    currentSampleData.GeneUnclear.add(data.Gene);
            }

            LOGGER.info("loaded {} driver gene records", mDriveGeneData.size());

        }
        catch (IOException exception)
        {
            LOGGER.error("Failed to read driver genes file({})", filename);
        }
    }

    private void loadSampleCountsData(final String filename)
    {
        if (filename.isEmpty())
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            // skip field names
            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("Empty sample counts CSV file({})", filename);
                return;
            }

            SampleData currentSampleData = null;

            while ((line = fileReader.readLine()) != null)
            {
                // parse CSV data
                String[] items = line.split(",");

                if (items.length != 5)
                    continue;

                SampleCountsData data = new SampleCountsData(
                        items[0], items[1], items[2], items[3], Integer.parseInt(items[4]));

                mSampleCountsData.add(data);

                if(!mCancerTypes.contains(data.CancerType))
                    mCancerTypes.add(data.CancerType);

                if(!mSamples.contains(data.SampleId))
                    mSamples.add(data.SampleId);

                if(!mCategories.contains(data.Category))
                {
                    mCategoryIndexMap.put(data.Category, mCategories.size());
                    mCategories.add(data.Category);
                }

                if(currentSampleData == null || !currentSampleData.SampleId.equals(data.SampleId))
                {
                    currentSampleData = getOrCreateSampleData(data.CancerType, data.SampleId);
                }

                if(data.Enriched.equals(VALUE_TRUE))
                    currentSampleData.CategoryKnown.add(data.Category);
                else if(data.Enriched.equals(VALUE_UNCLEAR))
                    currentSampleData.CategoryUnclear.add(data.Category);
            }

            LOGGER.info("loaded {} sample counts records", mSampleCountsData.size());

        }
        catch (IOException exception)
        {
            LOGGER.error("Failed to read sample counts file({})", filename);
        }
    }
}

class SampleData
{
    final String SampleId;

    List<String> GeneKnown;
    List<String> GeneUnclear;
    List<String> CategoryKnown;
    List<String> CategoryUnclear;

    public SampleData(final String sampleId)
    {
        SampleId = sampleId;
        GeneKnown = Lists.newArrayList();
        GeneUnclear = Lists.newArrayList();
        CategoryKnown = Lists.newArrayList();
        CategoryUnclear = Lists.newArrayList();
    }
}


class DriveGeneData
{
    public final String SampleId;
    public final String CancerType;
    public final String Gene;
    public final String DriverStatus;

    public DriveGeneData(final String sampleId, final String cancerType, final String gene, final String driverStatus)
    {
        SampleId = sampleId;
        CancerType = cancerType;
        Gene = gene;
        DriverStatus = driverStatus;
    }

}

class SampleCountsData
{
    public final String SampleId;
    public final String CancerType;
    public final String Category;
    public final String Enriched;
    public final int Count;

    public SampleCountsData(final String sampleId, final String cancerType, final String category, final String enriched, int count)
    {
        SampleId = sampleId;
        CancerType = cancerType;
        Category = category;
        Enriched = enriched;
        Count = count;
    }

}