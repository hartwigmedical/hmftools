package com.hartwig.hmftools.cup.common;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.VectorUtils.optimisedAdd;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.CuppaConfig.SUBSET_DELIM;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_NOISE_MEDIANS;
import static com.hartwig.hmftools.cup.common.ClassifierType.ALT_SJ_COHORT;
import static com.hartwig.hmftools.cup.common.ClassifierType.GENOMIC_POSITION_SIMILARITY;
import static com.hartwig.hmftools.cup.common.ClassifierType.SNV_96_PAIRWISE_SIMILARITY;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;

public class NoiseRefCache
{
    private final String mRefNoiseFile;
    private final Map<ClassifierType,double[]> mClassifierNoiseMedians;
    private final Map<ClassifierType,Integer> mClassifierNoiseAllocations;

    public static final String NOISE_ALLOC_DEFAULTS = "DEFAULTS";
    public static final String NOISE_ALLOC_NONE = "NONE";

    public NoiseRefCache(final String refDataDir)
    {
        mClassifierNoiseMedians = Maps.newHashMap();
        mClassifierNoiseAllocations = Maps.newHashMap();

        if(refDataDir != null)
        {
            mRefNoiseFile = refDataDir + REF_FILE_NOISE_MEDIANS;

            if(Files.exists(Paths.get(mRefNoiseFile)))
                loadAdjustments(mRefNoiseFile);
        }
        else
        {
            mRefNoiseFile = null;
        }
    }

    public boolean makeNoiseAdjustment(final ClassifierType classifierType)
    {
        return mClassifierNoiseMedians.containsKey(classifierType) && mClassifierNoiseAllocations.containsKey(classifierType);
    }

    public double[] getNoiseData(final ClassifierType classifierType) { return mClassifierNoiseMedians.get(classifierType); }
    public int getNoiseAllocation(final ClassifierType classifierType) { return mClassifierNoiseAllocations.get(classifierType); }

    public boolean hasNoiseData(final ClassifierType classifierType) { return mClassifierNoiseMedians.containsKey(classifierType); }
    public boolean hasNoiseAllocation(final ClassifierType classifierType) { return mClassifierNoiseAllocations.containsKey(classifierType); }

    public void addNoiseData(final ClassifierType classifierType, final double[] noiseMedians)
    {
        mClassifierNoiseMedians.put(classifierType, noiseMedians);
    }

    public void loadNoiseAllocations(final String config)
    {
        if(config == null)
            return;

        if(config.equals(NOISE_ALLOC_NONE))
            return;

        if(config.equals(NOISE_ALLOC_DEFAULTS))
        {
            mClassifierNoiseAllocations.put(ALT_SJ_COHORT, CupConstants.ALT_SJ_NOISE_ALLOCATION);
            mClassifierNoiseAllocations.put(SNV_96_PAIRWISE_SIMILARITY, CupConstants.SNV_96_NOISE_ALLOCATION);
            mClassifierNoiseAllocations.put(GENOMIC_POSITION_SIMILARITY, CupConstants.GEN_POS_COHORT_NOISE_ALLOCATION);

            // to be determined
            //mClassifierNoiseAllocations.put(GENOMIC_POSITION_PAIRWISE, CupConstants.GEN_POS_PAIRWISE_NOISE_ALLOCATION);
            return;
        }

        String[] classifierAllocations = config.split(SUBSET_DELIM, -1);

        for(String classifierStr : classifierAllocations)
        {
            String[] classifierValues = classifierStr.split("=", 2);
            mClassifierNoiseAllocations.put(ClassifierType.valueOf(classifierValues[0]), Integer.parseInt(classifierValues[1]));
        }
    }

    private void loadAdjustments(final String filename)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());
            fileData.remove(0); // ignore header

            for(final String line : fileData)
            {
                String[] values = line.split(DATA_DELIM, -1);
                ClassifierType classifierType = ClassifierType.valueOf(values[0]);

                double[] noiseAdjustments = new double[values.length - 1];

                for(int i = 1; i < values.length; ++i)
                {
                    noiseAdjustments[i - 1] = Double.parseDouble(values[i]);
                }

                mClassifierNoiseMedians.put(classifierType, noiseAdjustments);
            }

            if(!mClassifierNoiseMedians.isEmpty())
                CUP_LOGGER.info("loaded ref noise medians for classifiers: {}", mClassifierNoiseMedians.keySet());

        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read ref noise adjustments data file({}): {}", filename, e.toString());
        }
    }

    public void writeNoiseAdjustments()
    {
        if(mClassifierNoiseMedians.isEmpty() || mRefNoiseFile == null)
            return;

        try
        {
            BufferedWriter writer = createBufferedWriter(mRefNoiseFile, false);

            writer.write("Classifier,Medians");
            writer.newLine();

            for(Map.Entry<ClassifierType,double[]> entry : mClassifierNoiseMedians.entrySet())
            {
                writer.write(String.format("%s", entry.getKey()));

                final double[] adjustments = entry.getValue();

                for(int i = 0; i < adjustments.length; ++i)
                {
                    writer.write(String.format(",%.3f", adjustments[i]));
                }

                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write ref noise adjustments: {}", e.toString());
        }
    }

    public static double[] generateMedianValues(final Matrix matrix)
    {
        // determine the median count per bucket across the cancer types
        double[] bucketMedians = new double[matrix.Cols];

        final double[][] sourceData = matrix.getData();

        int dataItemCount = matrix.Rows; // ie samples, cancer types
        List<Double> sortedCounts = Lists.newArrayListWithCapacity(dataItemCount);

        int medianIndex = dataItemCount / 2;
        int medianLowerIndex = (dataItemCount % 2) == 0 ? medianIndex - 1 : medianIndex;

        if((dataItemCount % 2) == 1)
            ++medianIndex;

        for(int bucketIndex = 0; bucketIndex < matrix.Cols; ++bucketIndex)
        {
            sortedCounts.clear();

            for(int i = 0; i < matrix.Rows; ++i)
            {
                double count = sourceData[i][bucketIndex];
                optimisedAdd(sortedCounts, count, true);
            }

            double medianCount = 0;

            if(medianIndex == medianLowerIndex)
                medianCount = sortedCounts.get(medianIndex);
            else
                medianCount = (sortedCounts.get(medianLowerIndex) + sortedCounts.get(medianIndex)) * 0.5;

            bucketMedians[bucketIndex] = medianCount;
        }

        return bucketMedians;
    }

    public static void applyNoise(final Matrix matrix, final double[] bucketMedians, int noiseAllocation)
    {
        double medianTotal = sumVector(bucketMedians);

        // now scale these to the noise allocation
        final double[][] data = matrix.getData();
        for(int bucketIndex = 0; bucketIndex < matrix.Cols; ++bucketIndex)
        {
            double bucketMedian = bucketMedians[bucketIndex];
            double bucketPerc = bucketMedian / medianTotal;
            double bucketAlloc = bucketPerc * noiseAllocation;

            for(int i = 0; i < matrix.Rows; ++i)
            {
                data[i][bucketIndex] += bucketAlloc;
            }
        }
    }
}
