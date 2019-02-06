package com.hartwig.hmftools.svanalysis.stats;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.svanalysis.stats.GenericSampleData.SAMPLE_CAT_1_INDEX;
import static com.hartwig.hmftools.svanalysis.stats.GenericSampleData.SAMPLE_CAT_2_INDEX;

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
    private List<String> mSamples;

    // 2-way comparison, currently specific to driver genes and category data by cancer type
    private List<String> mGenes;
    private List<String> mCancerTypes;
    private List<String> mCategories;

    private Map<String, List<SampleData>> mCancerSampleData;
    private Map<String, Integer> mCategoryIndexMap;
    private Map<String, Integer> mGeneIndexMap;

    // generic data structure for 3-way co-occurrence
    private Map<String, List<GenericSampleData>> mGroupingSampleGenericData;
    private String mGroupingField;
    private String mCategory1;
    private String mCategory2;
    private List<String> mGroupingValues;
    private List<String> mCat1Values;
    private List<String> mCat2Values;

    private int[][][] mSampleCountsMatrix;

    private FisherExactTest mFisherET;

    private BufferedWriter mWriter;

    private static final Logger LOGGER = LogManager.getLogger(StatisticRoutines.class);

    public StatisticRoutines()
    {
        mCancerTypes = Lists.newArrayList();
        mSamples = Lists.newArrayList();
        mGenes = Lists.newArrayList();

        mCategories = Lists.newArrayList();
        mCategoryIndexMap = new HashMap();
        mGeneIndexMap = new HashMap();
        mCancerSampleData = new HashMap();
        mSampleCountsMatrix = null;

        mGroupingSampleGenericData = new HashMap();
        mGroupingValues = Lists.newArrayList();
        mCat1Values = Lists.newArrayList();
        mCat2Values = Lists.newArrayList();
        mGroupingField = "";
        mCategory1 = "";
        mCategory2 = "";

        mFisherET = new FisherExactTest();
        mWriter = null;
    }

    private static String DRIVER_GENES_FILE = "driver_genes_file";
    private static String SAMPLE_COUNTS_FILE = "sample_counts_file";
    private static String SAMPLE_GENERIC_FILE = "sample_generic_file";
    private static String OUTPUT_FILE = "stats_results_file";

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(DRIVER_GENES_FILE, true, "Drive genes file");
        options.addOption(SAMPLE_COUNTS_FILE, true, "Sample counts file");
        options.addOption(SAMPLE_GENERIC_FILE, true, "Sample data with 3 generic categories file");
        options.addOption(OUTPUT_FILE, true, "Results file");
    }

    public boolean loadConfig(final CommandLine cmd)
    {
        boolean valid = true;

        if(cmd.hasOption(DRIVER_GENES_FILE) && cmd.hasOption(SAMPLE_COUNTS_FILE))
        {
            loadDriverGeneData(cmd.getOptionValue(DRIVER_GENES_FILE));
            loadSampleCountsData(cmd.getOptionValue(SAMPLE_COUNTS_FILE));
            valid = initialiseTwoVariableOutput(cmd.getOptionValue(OUTPUT_FILE));
        }

        if(cmd.hasOption(SAMPLE_GENERIC_FILE))
        {
            loadSampleGenericData(cmd.getOptionValue(SAMPLE_GENERIC_FILE));
            valid = initialiseGenericThreeVariableOutput(cmd.getOptionValue(OUTPUT_FILE));
        }

        return valid;
    }

    public void runGenericThreeVariableStatisitics()
    {
        // for each of the group fields, calculate co-occurrence for each of the 2 categories
        for(final String groupingValue : mGroupingValues)
        {
            final List<GenericSampleData> sampleDataList = mGroupingSampleGenericData.get(groupingValue);

            if(sampleDataList == null || sampleDataList.isEmpty())
                continue;

            int sampleCount = sampleDataList.size();

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

                    for(final GenericSampleData sampleData : sampleDataList)
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
                                ++withCat1WithCat2;
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

                    double fisherProb = calcFisherExact(withCat1WithCat2, noCat1WithCat2, withCat1NoCat2, noCat1NoCat2, expectedVal);

                    //writeResultsData(cancerType, gene, category, sampleCount, withGeneTotal, withCatTotal, fisherProb,
                    //        expectedVal, withCatWithGene, noCatWithGene, withCatNoGene, noCatNoGene);

                }

            }

            /*
              categoryList = unique(tsgData$LohType)
  categoryCount = length(categoryList)

  totalSampleCount = unique(tsgSampleData$SampleId)

  for(geneName in geneList)
  {
    subResults = data.frame(matrix(ncol = 9, nrow = 0))
    colnames(subResults) = c("CatType", "CancerType", "SampleCount", "CatSC", "CancerSC", "WithCatWithCancerSC", "NoCatNoCancerSC", "ExpectedSC", "FisherET")

    geneData = tsgSampleData %>% filter(Gene==geneName)
    sampleCount = n_distinct(geneData$SampleId)

    print(paste("gene=", geneName, ", sampleCount=", sampleCount, sep=''))

    for(catName in categoryList)
    {
      # pan-cancer rates for each category (within this gene)
      scWithCat = nrow(geneData %>% filter(LohType==catName) %>% group_by(SampleId) %>% count())

      for(cancerType in cancerTypesList)
      {
        cancerData = tsgSampleData %>% filter(CancerType==cancerType&Gene==geneName)

        # scCancerType = nrow(tsgData %>% group_by(SampleId) %>% count())

        catCancerSummary = geneData %>% group_by(SampleId) %>% summarise(WithCat=sum(LohType==catName), WithCancer=sum(CancerType==cancerType))

        scNoCatNoCancer = nrow(catCancerSummary %>% filter(WithCat==0&WithCancer==0))
        scWithCatWithCancer = nrow(catCancerSummary %>% filter(WithCat>0&WithCancer>0))
        scWithCatNoCancer = nrow(catCancerSummary %>% filter(WithCat>0&WithCancer==0))
        scNoCatWithCancer = nrow(catCancerSummary %>% filter(WithCat==0&WithCancer>0))

        scWithCancer = nrow(catCancerSummary %>% filter(WithCancer>0))
        # scWithCat = nrow(catCancerSummary %>% filter(WithCat>0))

        expectedCount = round(scWithCat/sampleCount*scWithCancer,4)

        if(scWithCatWithCancer < 0 | scNoCatWithCancer < 0 | scWithCatNoCancer < 0 | scNoCatNoCancer < 0)
        {
          print(paste("INVALID catName=", catName, " cancer=", cancerType, " cancerSC=", scWithCancer, " catSC=", scWithCat, sep=''))

          print(paste("withCancer=", scWithCancer, " noCatWithCancer=", scNoCatWithCancer, " withCatWithCancer=", scWithCatWithCancer, sep=''))
          print(paste("noCancer=", scNoCancer, " noCatNoCancer=", scNoCatNoCancer, " withCatNoCancer=", scWithCatNoCancer, sep=''))
          return (allResults)
        }

        fishMatrix = rbind(c(scWithCatWithCancer,scNoCatWithCancer), c(scWithCatNoCancer,scNoCatNoCancer))




             */





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
            List<GenericSampleData> sampleDataList = null;

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

                String[] values = new String[GENERIC_DATA_CSV_COUNT];
                values[0] = sampleId;
                values[1] = groupingValue;
                values[2] = cat1Value;
                values[3] = cat2Value;

                boolean found = false;
                for(final GenericSampleData sampleData : sampleDataList)
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
                    GenericSampleData sampleData = new GenericSampleData(sampleId);
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

    private boolean initialiseGenericThreeVariableOutput(final String outputFileName)
    {
        try
        {
            mWriter = createBufferedWriter(outputFileName, false);

            mWriter.write(String.format("%s,%s,%s,SampleCount", mGroupingField, mCategory1, mCategory2));

            mWriter.write(String.format(",With%sCount,With%sCount,ExpectedCount,FETProb",
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

    private static String SPEC_CANCER = "";
    // private static String SPEC_CANCER = "Eye";
    private static String SPEC_GENE = "";
    // private static String SPEC_GENE = "OR11H1";
    private static String SPEC_CATEGORY = "";
    // private static String SPEC_CATEGORY = "DUP_LT_100";

    public void runTwoVariableStatistics()
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

    private boolean initialiseTwoVariableOutput(final String outputFileName)
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

            int recordCount = 0;
            SampleData currentSampleData = null;

            while ((line = fileReader.readLine()) != null)
            {
                // parse CSV data
                String[] items = line.split(",");

                if (items.length != 4)
                    continue;

                ++recordCount;

                final String sampleId = items[0];
                final String cancerType = items[1];
                final String gene = items[2];
                final String driverStatus = items[3];

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
            SampleData currentSampleData = null;

            while ((line = fileReader.readLine()) != null)
            {
                // parse CSV data
                String[] items = line.split(",");

                if (items.length != 5)
                    continue;

                ++recordCount;

                final String sampleId = items[0];
                final String cancerType = items[1];
                final String category = items[2];
                final String enriched = items[3];

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


class GenericSampleData
{
    public static int SAMPLE_CAT_1_INDEX = 0;
    public static int SAMPLE_CAT_2_INDEX = 1;

    public final String SampleId;

    private List<String[]> mCategoryData;

    List<String> GeneKnown;
    List<String> GeneUnclear;
    List<String> CategoryKnown;
    List<String> CategoryUnclear;

    public GenericSampleData(final String sampleId)
    {
        SampleId = sampleId;
        mCategoryData = Lists.newArrayList();
    }

    public void addCategoryData(final String cat1, final String cat2)
    {
        String[] values = {cat1, cat2};
        mCategoryData.add(values);
    }

    public final List<String[]> getCategoryData() { return mCategoryData; }
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

