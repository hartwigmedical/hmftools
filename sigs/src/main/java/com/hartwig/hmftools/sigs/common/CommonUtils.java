package com.hartwig.hmftools.sigs.common;

import static java.lang.Math.abs;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.sigs.SigUtils.calcAbsDiffs;
import static com.hartwig.hmftools.common.sigs.SigUtils.calcLinearLeastSquares;
import static com.hartwig.hmftools.common.sigs.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.common.utils.MatrixUtils;

import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class CommonUtils
{
    public static final String SAMPLE_IDS = "sample";
    public static final String SAMPLE_COUNTS_FILE = "sample_counts_file";
    public static final String SIGNATURES_FILE = "signatures_file";
    public static final String LOG_DEBUG = "log_debug";

    public static final String OUTPUT_FILE_ID = "output_file_id";

    public static final Logger SIG_LOGGER = LogManager.getLogger(CommonUtils.class);

    public static void addCmdLineArgs(final Options options)
    {
        options.addOption(SAMPLE_IDS, true, "Optional - list of sampleIds, separated by ';");
        options.addOption(SAMPLE_COUNTS_FILE, true, "Path to the main input file");
        options.addOption(SIGNATURES_FILE, true, "Signature definitions");
        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(OUTPUT_FILE_ID, true, "Output file ID");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
    }

    public static String formOutputFilename(final String outputDir, final String outputId, final String fileId)
    {
        if(outputId != null)
            return outputDir + fileId + "." + outputId + ".csv";
        else
            return outputDir + fileId  + ".csv";
    }

    public static BufferedWriter getNewFile(final String outputDir, final String fileName) throws IOException
    {
        if(outputDir.isEmpty())
            return null;

        return createBufferedWriter(outputDir + fileName, false);
    }

    public static Matrix loadSampleMatrixCounts(final String filename, final List<String> sampleIds)
    {
        return MatrixUtils.loadMatrixDataFile(filename, sampleIds, Lists.newArrayList("BucketName"));
    }

    public static List<Integer> getCombinedList(final List<Integer> list1, final List<Integer> list2)
    {
        // gets super set, including non-common values
        List<Integer> combinedSet = Lists.newArrayList();

        for(Integer value : list1)
        {
            if(!combinedSet.contains(value))
                combinedSet.add(value);
        }

        for(Integer value : list2)
        {
            if(!combinedSet.contains(value))
                combinedSet.add(value);
        }

        return combinedSet;
    }

    public static List<Integer> getMatchingList(final List<Integer> list1, final List<Integer> list2)
    {
        // gets union/common set
        List<Integer> matchedList = Lists.newArrayList();

        for(Integer value : list1)
        {
            if(list2.contains(value))
                matchedList.add(value);
        }

        return matchedList;
    }

    public static List<Integer> getDiffList(final List<Integer> list1, final List<Integer> list2)
    {
        // returns list of values in 1 but not in 2
        List<Integer> diffList = Lists.newArrayList();

        for(Integer value : list1)
        {
            if(!list2.contains(value))
                diffList.add(value);
        }

        return diffList;
    }

    public static double calcMinPositiveRatio(final double[] params, final double[] data)
    {
        if(data.length != params.length)
            return 0;

        // returns the max ratio applying the params to the data
        // where the fit values does not exceed the actual data
        double minRatio = 0;

        for(int i = 0; i < data.length; ++i)
        {
            if(data[i] == 0)
                continue;

            double ratio = data[i] / params[i];

            if(ratio < minRatio || minRatio == 0)
                minRatio = ratio;
        }

        return minRatio;
    }

    public static double calcBestFitWithinProbability(int itemId, final double[] ratios, final double[] data,
            double requiredProb, double reqResidualsPerc)
    {
        if(data.length != ratios.length)
            return 0;

        int itemCount = data.length;

        // calculate least squares and min positive as the upper and lower starting bounds
        double lsAlloc = calcLinearLeastSquares(ratios, data);
        double minPosAlloc = calcMinPositiveRatio(ratios, data);
        double dataTotal = sumVector(data);

        double[] currentFit = new double[data.length];
        double[] reducedData = new double[data.length];

        int iterations = 0;
        int maxIterations = 10;

        double currentAlloc = dataTotal;
        double lowerAlloc = minPosAlloc;
        double upperAlloc = dataTotal;
        double probDiff = 0;
        double currentProb = 0;
        double residuals = 0;
        double residualsPerc = 0;

        while(iterations < maxIterations)
        {
            // work out probability of the current fit
            for(int i = 0; i < itemCount; ++i)
            {
                currentFit[i] = currentAlloc * ratios[i];
                reducedData[i] = min(currentFit[i], data[i]);
            }

            residuals = calcAbsDiffs(currentFit, reducedData);
            currentProb = CssRoutines.calcLogLikelihood(currentFit, reducedData, false);
            probDiff = abs(requiredProb - currentProb) / requiredProb;
            residualsPerc = residuals / dataTotal;

            if(probDiff < 0.1 && residualsPerc < reqResidualsPerc)
                break;

            // if prob is too high, need to increase the counts allocation
            if(currentProb > requiredProb)
            {
                if(currentAlloc >= upperAlloc - 1)
                    break;

                lowerAlloc = currentAlloc;
                currentAlloc = (int)round((currentAlloc + upperAlloc) * 0.5);
            }
            else
            {
                if(currentAlloc <= lowerAlloc + 1)
                    break;

                upperAlloc = currentAlloc;
                currentAlloc = (int)round((currentAlloc + lowerAlloc) * 0.5);
            }

            ++iterations;
        }

        SIG_LOGGER.debug(String.format("sample(%d) total(%.0f) finalAlloc(%.0f minRatio=%.0f leastSq=%.0f) residuals(%.0f perc=%.3f) prob(%.4f) iter(%s)",
                itemId, dataTotal, currentAlloc, minPosAlloc, lsAlloc, residuals, residualsPerc, currentProb,
                iterations >= maxIterations ? "max" : String.valueOf(iterations)));

        return currentAlloc;
    }

    public static List<String> loadSampleListFile(final String filename)
    {
        try
        {
            final List<String> sampleIds = Files.readAllLines(new File(filename).toPath());

            if (sampleIds.get(0).equals("SampleId"))
                sampleIds.remove(0);

            SIG_LOGGER.info("Loaded {} specific sample IDs", sampleIds.size());

            return sampleIds;
        }
        catch (IOException exception)
        {
            SIG_LOGGER.error("failed to read sample list input CSV file({}): {}", filename, exception.toString());
            return Lists.newArrayList();
        }
    }

}
