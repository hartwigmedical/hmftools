package com.hartwig.hmftools.isofox.expression.cohort;

import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.stats.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.stats.Percentiles.getPercentile;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

public class CohortGenePercentiles
{
    private final Map<String,Map<String,double[]>> mGenePercentiles; // by geneId and then cancer type
    private final Map<String,Map<String,Double>> mGeneMedians;

    public static final String PAN_CANCER = "ALL";
    public static final String CANCER_TYPE_OTHER = "Other";

    public CohortGenePercentiles(final String cohortFile)
    {
        mGenePercentiles = Maps.newHashMap();
        mGeneMedians = Maps.newHashMap();

        loadCohortFile(cohortFile);
    }

    private void loadCohortFile(final String filename)
    {
        if(filename == null || !Files.exists(Paths.get(filename)))
        {
            ISF_LOGGER.error("invalid gene cohort distribution file({})", filename);
            return;
        }

        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            int geneIdIndex = fieldsIndexMap.get(FLD_GENE_ID);
            int cancerIndex = fieldsIndexMap.get("CancerType");
            int medianIndex = fieldsIndexMap.get("Median");
            int percStartIndex = fieldsIndexMap.get("Pct_0");

            String currentGeneId = "";
            Map<String,double[]> percentilesMap = null;
            Map<String,Double> medianMap = null;

            for(String line : lines)
            {
                final String[] items = line.split(DELIMITER, -1);

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

            ISF_LOGGER.info("loaded distribution for {} genes from file({})", mGeneMedians.size(), filename);
        }
        catch (IOException e)
        {
            ISF_LOGGER.warn("failed to load gene distribution file({}): {}", filename, e.toString());
        }
    }

    public double getTpmMedian(final String geneId, final String cancerType)
    {
        final Map<String,Double> medianMap = mGeneMedians.get(geneId);

        if(medianMap == null || !medianMap.containsKey(cancerType))
            return -1;

        return medianMap.get(cancerType);
    }

    public double getTpmPercentile(final String geneId, final String cancerType, final double sampleTpm)
    {
        final Map<String,double[]> percentileMap = mGenePercentiles.get(geneId);

        if(percentileMap == null || !percentileMap.containsKey(cancerType))
            return -1;

        return getPercentile(percentileMap.get(cancerType), sampleTpm);
    }
}
