package com.hartwig.hmftools.statcalcs.cooc;

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

/* 2-way comparison, currently specific to driver genes and category data by cancer type

    input files:
        1. SampleCounts - SampleId,CancerType,Category,Enriched,Count
        2. DriverGenes - SampleId,CancerType,Gene,DriverStatus
*/


public class SampleCountsCoOccurence
{
    private List<String> mSamples;
    private Map<String, List<SampleGeneData>> mCancerSampleData;

    private List<String> mGenes;
    private List<String> mCancerTypes;
    private List<String> mCategories;

    private Map<String, Integer> mCategoryIndexMap;
    private Map<String, Integer> mGeneIndexMap;

    private int[][][] mSampleCountsMatrix;

    private final FisherExactTest mFisherET;

    private BufferedWriter mWriter;

    private static final String DRIVER_GENES_FILE = "driver_genes_file";
    private static final String SAMPLE_COUNTS_FILE = "sample_counts_file";

    private static final Logger LOGGER = LogManager.getLogger(SampleCountsCoOccurence.class);

    public SampleCountsCoOccurence(final CommandLine cmd, final String outputDir)
    {
        mCancerTypes = Lists.newArrayList();
        mSamples = Lists.newArrayList();
        mGenes = Lists.newArrayList();

        mCategories = Lists.newArrayList();
        mCategoryIndexMap = new HashMap<>();
        mGeneIndexMap = new HashMap<>();
        mCancerSampleData = new HashMap<>();
        mSampleCountsMatrix = null;

        mFisherET = new FisherExactTest();

        final String sampleCountsFile = cmd.getOptionValue(SAMPLE_COUNTS_FILE);
        final String driverDataFile = cmd.getOptionValue(DRIVER_GENES_FILE);

        loadDriverGeneData(driverDataFile);
        loadSampleCountsData(sampleCountsFile);

        final String outputFile = outputDir + "STATS_GENE_CATEGORY.csv";
        initialiseOutputFile(outputFile);
    }

    public static void addCmdLineOptions(Options options)
    {
        options.addOption(DRIVER_GENES_FILE, true, "Drive genes file");
        options.addOption(SAMPLE_COUNTS_FILE, true, "Sample counts file");
    }

    public static boolean hasConfig(final CommandLine cmd)
    {
        return cmd.hasOption(DRIVER_GENES_FILE) && cmd.hasOption(SAMPLE_COUNTS_FILE);
    }

    // non-generic 2-variable test, using 2 distinct data sets

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

    // debug
    private static String SPEC_CANCER = "";
    // private static String SPEC_CANCER = "Eye";
    private static String SPEC_GENE = "";
    // private static String SPEC_GENE = "OR11H1";
    private static String SPEC_CATEGORY = "";
    // private static String SPEC_CATEGORY = "DUP_LT_100";

    public void run()
    {
        if(mCancerSampleData.isEmpty())
            return;

        List<SampleGeneData> allSampleDataList = Lists.newArrayList();

        for (final String cancerType : mCancerTypes)
        {
            if (!SPEC_CANCER.isEmpty() && !cancerType.equals(SPEC_CANCER))
                continue;

            List<SampleGeneData> sampleDataList = mCancerSampleData.get(cancerType);
            allSampleDataList.addAll(sampleDataList);

            analyseCancerType(cancerType, sampleDataList);

        }

        // repeat for all cancers combined
        analyseCancerType("All", allSampleDataList);

        FileWriterUtils.closeBufferedWriter(mWriter);

        LOGGER.info("analysis complete");
    }

    private void analyseCancerType(final String cancerType, final List<SampleGeneData> sampleDataList)
    {
        /* calc method (for each cancer type):
        - each sample data has enriched categories and driver genes
        - from this it can contribute to each non-false cell
        - also tally up across all samples the counts per category and gene
         */
        int categoryCount = mCategories.size();
        int geneCount = mGenes.size();

        int hypothesesCount = geneCount * categoryCount;

        int sampleCount = sampleDataList.size();
        LOGGER.info("processing cancerType({}) with {} samples", cancerType, sampleCount);

        // populate the matrix with counts
        mSampleCountsMatrix = new int[geneCount][categoryCount][INDEX_COUNT];

        mFisherET.initialise(sampleCount);

        Map<String, Integer> withCategoryTotals = new HashMap<>();
        Map<String, Integer> unclearCategoryTotals = new HashMap<>();
        Map<String, Integer> withGeneTotals = new HashMap<>();
        Map<String, Integer> unclearGeneTotals = new HashMap<>();

        for(final SampleGeneData sampleData : sampleDataList)
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

                double geneSamplesPerc = withGeneTotal/(double)sampleCount;
                double expectedVal  = withCatTotal * geneSamplesPerc;
                double fisherProb = mFisherET.calc(withCatWithGene, noCatWithGene, withCatNoGene, noCatNoGene, expectedVal);

                writeResultsData(cancerType, gene, category, sampleCount, withGeneTotal, withCatTotal, fisherProb,
                        expectedVal, hypothesesCount, withCatWithGene, noCatWithGene, withCatNoGene, noCatNoGene);
            }
        }

        LOGGER.info("cancerType({}) results written to file", cancerType);
    }

    private boolean initialiseOutputFile(final String outputFileName)
    {
        try
        {
            mWriter = FileWriterUtils.createBufferedWriter(outputFileName, false);

            mWriter.write("CancerType,Gene,Category,SampleCount");
            mWriter.write(",WithGeneCount,WithCategoryCount,ExpectedCount,FETProb");
            mWriter.write(",WithCatWithGene,NoCatWithGene,WithCatNoGene,NoCatNoGene,TestCount,CountGtExp");

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
            int withGeneTotal, int withCatTotal, double fetProbability, double expectedVal, int hypothesesCount,
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
                    String.format(",%d,%d,%.2f,%4.3e,%d,%d,%d,%d,%d,%s",
                            withGeneTotal, withCatTotal, expectedVal, fetProbability,
                            withCatWithGene, noCatWithGene, withCatNoGene, noCatNoGene,
                            hypothesesCount, withCatWithGene > expectedVal));

            mWriter.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to stats output file: {}", e.toString());
        }
    }
    private SampleGeneData getOrCreateSampleData(final String cancerType, final String sampleId)
    {
        List<SampleGeneData> sampleDataList = mCancerSampleData.get(cancerType);

        SampleGeneData sampleData = null;
        if (sampleDataList == null)
        {
            sampleDataList = Lists.newArrayList();
            mCancerSampleData.put(cancerType, sampleDataList);
        }
        else
        {
            for (final SampleGeneData sd : sampleDataList)
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
            sampleData = new SampleGeneData(sampleId);
            sampleDataList.add(sampleData);
        }

        return sampleData;
    }

    private static final int SC_COL_SAMPLEID = 0;
    private static final int SC_COL_CT = 1;
    private static final int SC_COL_GENE = 2;
    private static final int SC_COL_DS = 3;

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

            int recordCount = 0;
            SampleGeneData currentSampleData = null;

            while ((line = fileReader.readLine()) != null)
            {
                // parse CSV data
                String[] items = line.split(",");

                if (items.length != 4)
                    continue;

                ++recordCount;

                final String sampleId = items[SC_COL_SAMPLEID];
                final String cancerType = items[SC_COL_CT];
                final String gene = items[SC_COL_GENE];
                final String driverStatus = items[SC_COL_DS];

                if(!mCancerTypes.contains(cancerType))
                    mCancerTypes.add(cancerType);

                if(!mSamples.contains(sampleId))
                    mSamples.add(sampleId);

                if(!mGenes.contains(gene))
                {
                    mGeneIndexMap.put(gene, mGenes.size());
                    mGenes.add(gene);
                }

                if(currentSampleData == null || !currentSampleData.SampleId.equals(sampleId))
                {
                    currentSampleData = getOrCreateSampleData(cancerType, sampleId);
                }

                if(driverStatus.equals(VALUE_TRUE))
                    currentSampleData.GeneKnown.add(gene);
                else if(driverStatus.equals(VALUE_UNCLEAR))
                    currentSampleData.GeneUnclear.add(gene);
            }

            LOGGER.info("loaded {} driver gene records", recordCount);

        }
        catch (IOException exception)
        {
            LOGGER.error("Failed to read driver genes file({})", filename);
        }
    }

    private static final int SC_COL_CAT = 2;
    private static final int SC_COL_ENRICHED = 3;

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

            int recordCount = 0;
            SampleGeneData currentSampleData = null;

            while ((line = fileReader.readLine()) != null)
            {
                // parse CSV data
                String[] items = line.split(",");

                if (items.length != 5)
                    continue;

                ++recordCount;

                final String sampleId = items[SC_COL_SAMPLEID];
                final String cancerType = items[SC_COL_CT];
                final String category = items[SC_COL_CAT];
                final String enriched = items[SC_COL_ENRICHED];

                if(!mCancerTypes.contains(cancerType))
                    mCancerTypes.add(cancerType);

                if(!mSamples.contains(sampleId))
                    mSamples.add(sampleId);

                if(!mCategories.contains(category))
                {
                    mCategoryIndexMap.put(category, mCategories.size());
                    mCategories.add(category);
                }

                if(currentSampleData == null || !currentSampleData.SampleId.equals(sampleId))
                {
                    currentSampleData = getOrCreateSampleData(cancerType, sampleId);
                }

                if(enriched.equals(VALUE_TRUE))
                    currentSampleData.CategoryKnown.add(category);
                else if(enriched.equals(VALUE_UNCLEAR))
                    currentSampleData.CategoryUnclear.add(category);
            }

            LOGGER.info("loaded {} sample counts records", recordCount);

        }
        catch (IOException exception)
        {
            LOGGER.error("Failed to read sample counts file({})", filename);
        }
    }


}
