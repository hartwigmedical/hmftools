package com.hartwig.hmftools.common.isofox;

import static com.hartwig.hmftools.common.stats.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.stats.Percentiles.getPercentile;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class GeneExpressionDistributionData
{
    private static final Logger LOGGER = LogManager.getLogger(GeneExpressionDistributionData.class);
    public static final int NOT_AVAILABLE = -1;

    private final Map<String, Map<String, double[]>> mGenePercentiles; // by geneId and then cancer type
    private final Map<String, Map<String, Double>> mGeneMedians;

    public static final String PAN_CANCER = "ALL";

    public GeneExpressionDistributionData(final String cohortFile) throws IOException
    {
        mGenePercentiles = Maps.newHashMap();
        mGeneMedians = Maps.newHashMap();

        loadCohortFile(cohortFile);
    }

    private void loadCohortFile(final String filename) throws IOException
    {
        final List<String> lines = Files.readAllLines(Paths.get(filename));

        String fileDelim = inferFileDelimiter(filename);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), fileDelim);
        lines.remove(0);

        int geneIdIndex = fieldsIndexMap.get(FLD_GENE_ID);
        int cancerIndex = fieldsIndexMap.get("CancerType");
        int medianIndex = fieldsIndexMap.get("Median");
        int percStartIndex = fieldsIndexMap.get("Pct_0");

        String currentGeneId = "";
        Map<String, double[]> percentilesMap = null;
        Map<String, Double> medianMap = null;

        for(String line : lines)
        {
            final String[] items = line.split(fileDelim, -1);

            final String geneId = items[geneIdIndex];
            final String cancerType = items[cancerIndex];
            final double median = Double.parseDouble(items[medianIndex]);

            if(!currentGeneId.equals(geneId))
            {
                currentGeneId = geneId;
                percentilesMap = Maps.newHashMap();
                medianMap = Maps.newHashMap();
                mGeneMedians.put(geneId, medianMap);
                mGenePercentiles.put(geneId, percentilesMap);
            }

            medianMap.put(cancerType, median);

            double[] percentileData = new double[PERCENTILE_COUNT];

            for(int i = percStartIndex; i < items.length; ++i)
            {
                double tpm = Double.parseDouble(items[i]);
                percentileData[i - percStartIndex] = tpm;
            }

            percentilesMap.put(cancerType, percentileData);
        }

        LOGGER.info(" Loaded distribution for {} genes from {}", mGeneMedians.size(), filename);
    }

    public double getTpmMedian(final String geneId, final String cancerType)
    {
        final Map<String, Double> medianMap = mGeneMedians.get(geneId);

        if(medianMap == null || !medianMap.containsKey(cancerType))
        {
            return NOT_AVAILABLE;
        }

        return medianMap.get(cancerType);
    }

    public double getTpmPercentile(final String geneId, final String cancerType, final double sampleTpm)
    {
        final Map<String, double[]> percentileMap = mGenePercentiles.get(geneId);

        if(percentileMap == null || !percentileMap.containsKey(cancerType))
        {
            return NOT_AVAILABLE;
        }

        return getPercentile(percentileMap.get(cancerType), sampleTpm);
    }

    Set<String> configuredCancerTypes()
    {
        Set<String> configuredCancerTypes = Sets.newHashSet();
        for(Map<String, Double> geneMedian : mGeneMedians.values())
        {
            configuredCancerTypes.addAll(geneMedian.keySet());
        }
        return configuredCancerTypes;
    }
}
