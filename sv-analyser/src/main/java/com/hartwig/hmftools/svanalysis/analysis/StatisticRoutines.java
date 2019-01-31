package com.hartwig.hmftools.svanalysis.analysis;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svanalysis.annotators.DriverGeneAnnotator;

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

    private Map<String, List<String>> mCancerSamples;

    private int[][][] mSampleCountsMatrix;

    private static final Logger LOGGER = LogManager.getLogger(StatisticRoutines.class);

    public StatisticRoutines()
    {
        mDriveGeneData = Lists.newArrayList();
        mSampleCountsData = Lists.newArrayList();
        mCancerTypes = Lists.newArrayList();
        mSamples = Lists.newArrayList();
        mGenes = Lists.newArrayList();
        mCategories = Lists.newArrayList();
        mCancerSamples = new HashMap();
        mSampleCountsMatrix = null;
    }

    private static String DRIVER_GENES_FILE = "driver_genes_file";
    private static String SAMPLE_COUNTS_FILE = "sample_counts_file";

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(DRIVER_GENES_FILE, true, "Drive genes file");
        options.addOption(SAMPLE_COUNTS_FILE, true, "Sample counts file");
    }

    public boolean loadConfig(final CommandLine cmd)
    {
        loadDriverGeneData(cmd.getOptionValue(DRIVER_GENES_FILE));
        loadSampleCountsData(cmd.getOptionValue(SAMPLE_COUNTS_FILE));

        return true;
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

    public void runStatistics()
    {
        /*
        int categoryCount = mCategories.size();
        int geneCount = mGenes.size();
        mSampleCountsMatrix = new int[geneCount][categoryCount][INDEX_COUNT];

        for(final String cancerType : mCancerTypes)
        {
            List<String> samples = mCancerSamples.get(cancerType);

            int sampleCount = samples.size();
            LOGGER.debug("processing cancerType({}) with {} samples", cancerType, sampleCount);

            // find
            for(final SampleCountsData countsData : mSampleCountsData)
            {
                if (!countsData.CancerType.equals(cancerType) || !countsData.Category.equals(category))
                    continue;

                for (final DriveGeneData dgData : mDriveGeneData)
                {
                    if (!dgData.CancerType.equals(cancerType) || !dgData.Gene.equals(gene))
                        continue;

                    if (!countsData.SampleId.equals(dgData.SampleId))
                        continue;

                    if (countsData.Enriched.equals(VALUE_TRUE))
                    {
                        if (dgData.DriverStatus.equals(VALUE_TRUE))
                            ++scCatEnrichedGeneDriver;
                        else if (dgData.DriverStatus.equals(VALUE_UNCLEAR))
                            ++scCatEnrichedGeneUnclear;
                    }
                    else if (countsData.Enriched.equals(VALUE_UNCLEAR))
                    {
                        if (dgData.DriverStatus.equals(VALUE_TRUE))
                            ++scCatUnclearGeneDriver;
                        else if (dgData.DriverStatus.equals(VALUE_UNCLEAR))
                            ++scCatUnclearGeneUnclear;
                    }
                }
            }


                    for(final String category : mCategories)
            {
                // count up samples with this category
                int scCatEnriched = 0;
                int scCatUnclear = 0;

                for(final SampleCountsData countsData : mSampleCountsData)
                {
                    if(!countsData.CancerType.equals(cancerType) || !countsData.Category.equals(category))
                        continue;

                    if(countsData.Enriched.equals(VALUE_TRUE))
                        ++scCatEnriched;
                    else if(countsData.Enriched.equals(VALUE_UNCLEAR))
                        ++scCatUnclear;
                }

                int scCatNotEnriched = sampleCount - scCatEnriched - scCatUnclear;

                for(final String gene : mGenes)
                {
                    int scGeneDriver = 0;
                    int scGeneUnclear = 0;

                    for(final DriveGeneData dgData : mDriveGeneData)
                    {
                        if(!dgData.CancerType.equals(cancerType) || !dgData.Gene.equals(gene))
                            continue;

                        if(dgData.DriverStatus.equals(VALUE_TRUE))
                            ++scGeneDriver;
                        else if(dgData.DriverStatus.equals(VALUE_UNCLEAR))
                            ++scGeneUnclear;
                    }

                    int scGeneNotDriver = sampleCount - scGeneDriver - scGeneUnclear;

                    // work out how frequently each sample appears in each gene and category
                    int scCatEnrichedGeneDriver = 0;
                    int scCatEnrichedGeneUnclear = 0;
                    int scCatUnclearGeneDriver = 0;
                    int scCatUnclearGeneUnclear = 0;

                    for(final SampleCountsData countsData : mSampleCountsData)
                    {
                        if (!countsData.CancerType.equals(cancerType) || !countsData.Category.equals(category))
                            continue;

                        for (final DriveGeneData dgData : mDriveGeneData)
                        {
                            if (!dgData.CancerType.equals(cancerType) || !dgData.Gene.equals(gene))
                                continue;

                            if(!countsData.SampleId.equals(dgData.SampleId))
                                continue;

                            if(countsData.Enriched.equals(VALUE_TRUE))
                            {
                                if(dgData.DriverStatus.equals(VALUE_TRUE))
                                    ++scCatEnrichedGeneDriver;
                                else if(dgData.DriverStatus.equals(VALUE_UNCLEAR))
                                    ++scCatEnrichedGeneUnclear;
                            }
                            else if(countsData.Enriched.equals(VALUE_UNCLEAR))
                            {
                                if(dgData.DriverStatus.equals(VALUE_TRUE))
                                    ++scCatUnclearGeneDriver;
                                else if(dgData.DriverStatus.equals(VALUE_UNCLEAR))
                                    ++scCatUnclearGeneUnclear;
                            }

                            break;
                        }
                    }

                    // now calculate

                }


            }


        }
        */
    }

    private void addCancerSample(final String cancerType, final String sampleId)
    {
        List<String> sampleList = mCancerSamples.get(cancerType);

        if(sampleList == null)
        {
            sampleList = Lists.newArrayList();
            mCancerSamples.put(cancerType, sampleList);
        }

        if(!sampleList.contains(sampleId))
            sampleList.add(sampleId);
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
                    mGenes.add(data.Gene);

                addCancerSample(data.CancerType, data.SampleId);
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
                    mCategories.add(data.Category);

                addCancerSample(data.CancerType, data.SampleId);
            }

            LOGGER.info("loaded {} sample counts records", mSampleCountsData.size());

        }
        catch (IOException exception)
        {
            LOGGER.error("Failed to read sample counts file({})", filename);
        }
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