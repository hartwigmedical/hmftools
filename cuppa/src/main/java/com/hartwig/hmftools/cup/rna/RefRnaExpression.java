package com.hartwig.hmftools.cup.rna;

import static java.lang.Math.log;

import static com.hartwig.hmftools.common.sigs.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.ref.RefDataConfig;

public class RefRnaExpression
{
    private final RefDataConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private final Map<String,Map<String,Double>> mGeneCancerExpressionData;
    private final Map<String,String> mGeneIdNameMap;

    public RefRnaExpression(final RefDataConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mGeneCancerExpressionData = Maps.newHashMap();
        mGeneIdNameMap = Maps.newHashMap();
    }

    public void buildRefDataSets()
    {
        if(mConfig.RefRnaGeneExpFile.isEmpty())
            return;

        loadRefRnaGeneExpression(mConfig.RefRnaGeneExpFile);

        try
        {
            final String filename = mConfig.OutputDir + "cup_ref_ct_rna_gene_exp.csv";
            BufferedWriter writer = createBufferedWriter(filename, false);

            final Set<String> cancerTypes = mSampleDataCache.RefCancerSampleData.keySet();

            writer.write("GeneId,GeneName");

            for(final String cancerType : cancerTypes)
            {
                writer.write(String.format(",%s", cancerType));
            }

            writer.newLine();

            for(Map.Entry<String,Map<String,Double>> geneEntry : mGeneCancerExpressionData.entrySet())
            {
                final String geneId = geneEntry.getKey();
                final String geneName = mGeneIdNameMap.get(geneId);

                writer.write(String.format("%s,%s", geneId, geneName));

                final Map<String,Double> cancerTpmTotals = geneEntry.getValue();

                for(final String cancerType : cancerTypes)
                {
                    Double tpmTotal = cancerTpmTotals.get(cancerType);

                    writer.write(String.format(",%.1f", tpmTotal != null ? tpmTotal : 0));
                }

                writer.newLine();
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write ref RNA gene expression output: {}", e.toString());
        }
    }

    private void loadRefRnaGeneExpression(final String filename)
    {
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            // SampleId,CancerType,GeneId,GeneName,TPM,CohortMedianTPM,CancerMedianTPM,CohortPercentile,CancerPercentile

            // SampleId,Chromosome,Position,Count
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(line, ",");

            int sampleIndex = fieldsIndexMap.get("SampleId");
            int geneIdIndex = fieldsIndexMap.get("GeneId");
            int geneNameIndex = fieldsIndexMap.get("GeneName");
            int tpmIndex = fieldsIndexMap.get("TPM");

            line = fileReader.readLine(); // skip header

            while (line != null)
            {
                final String[] items = line.split(",", -1);
                String sampleId = items[sampleIndex];

                final String cancerType = mSampleDataCache.RefSampleCancerTypeMap.get(sampleId);

                if(cancerType == null)
                {
                    line = fileReader.readLine();
                    continue;
                }

                String geneId = items[geneIdIndex];
                String geneName = items[geneNameIndex];
                double tpm = Double.parseDouble(items[tpmIndex]);

                Map<String,Double> cancerTpmTotals = mGeneCancerExpressionData.get(geneId);

                if(cancerTpmTotals == null)
                {
                    cancerTpmTotals = Maps.newHashMap();
                    mGeneCancerExpressionData.put(geneId, cancerTpmTotals);
                    mGeneIdNameMap.put(geneId, geneName);
                }

                double scaledTpm = scaleTpm(tpm);
                Double tpmTotal = cancerTpmTotals.get(cancerType);

                if(tpmTotal == null)
                    cancerTpmTotals.put(cancerType, scaledTpm);
                else
                    cancerTpmTotals.put(cancerType, tpmTotal + scaledTpm);

                line = fileReader.readLine();
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read RNA gene expression data file({}): {}", filename, e.toString());
        }
    }

    private double scaleTpm(double tpm)
    {
        if(tpm <= 1)
            return 0;

        return log(tpm);
    }
}
