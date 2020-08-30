package com.hartwig.hmftools.cup.rna;

import static java.lang.Math.log;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.ref.RefDataConfig;

public class RefRnaExpression
{
    private final RefDataConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private final Map<String,Map<String,Double>> mGeneCancerExpressionData;
    private final Map<String,Map<String,Double>> mGeneSampleExpressionData;
    private final Map<String,String> mGeneIdNameMap;
    private final List<String> mSampleIds;

    public RefRnaExpression(final RefDataConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mGeneCancerExpressionData = Maps.newHashMap();
        mGeneSampleExpressionData = Maps.newHashMap();
        mGeneIdNameMap = Maps.newHashMap();
        mSampleIds = Lists.newArrayList();
    }

    public void buildRefDataSets()
    {
        if(mConfig.RefRnaGeneExpFile.isEmpty())
            return;

        CUP_LOGGER.debug("loading RNA gene expression data");

        loadRefRnaGeneExpression(mConfig.RefRnaGeneExpFile);

        CUP_LOGGER.debug("writing RNA gene expression reference data");

        writeGeneData();
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

            int samplesProcessed = 0;
            String currentSample = "";

            while (line != null)
            {
                final String[] items = line.split(",", -1);
                String sampleId = items[sampleIndex];

                if(!currentSample.equals(sampleId))
                {
                    currentSample = sampleId;
                    mSampleIds.add(sampleId);
                    ++samplesProcessed;

                    if(samplesProcessed > 0 && (samplesProcessed % 100) == 0)
                    {
                        CUP_LOGGER.info("processed {} RNA samples", samplesProcessed);
                    }
                }

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

                // stored by sampleId
                Map<String,Double> sampleTpms = mGeneSampleExpressionData.get(geneId);

                if(sampleTpms == null)
                {
                    sampleTpms = Maps.newHashMap();
                    mGeneSampleExpressionData.put(geneId, sampleTpms);
                }

                sampleTpms.put(sampleId, scaledTpm);

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

    private void writeGeneData()
    {
        try
        {
            final String filename = mConfig.OutputDir + "cup_ref_ct_rna_gene_exp.csv";
            BufferedWriter cancerWriter = createBufferedWriter(filename, false);

            final String sampleDataFilename = mConfig.OutputDir + "cup_ref_rna_gene_exp.csv";
            BufferedWriter sampleWriter = createBufferedWriter(sampleDataFilename, false);

            final Set<String> cancerTypes = mSampleDataCache.RefCancerSampleData.keySet();

            cancerWriter.write("GeneId,GeneName");

            for(final String cancerType : cancerTypes)
            {
                cancerWriter.write(String.format(",%s", cancerType));
            }

            cancerWriter.newLine();

            sampleWriter.write("GeneId,GeneName");

            for(final String sampleId : mSampleIds)
            {
                sampleWriter.write(String.format(",%s", sampleId));
            }

            sampleWriter.newLine();

            for(Map.Entry<String,Map<String,Double>> geneEntry : mGeneCancerExpressionData.entrySet())
            {
                final String geneId = geneEntry.getKey();
                final String geneName = mGeneIdNameMap.get(geneId);

                cancerWriter.write(String.format("%s,%s", geneId, geneName));

                final Map<String,Double> cancerTpmTotals = geneEntry.getValue();

                for(final String cancerType : cancerTypes)
                {
                    Double tpmTotal = cancerTpmTotals.get(cancerType);
                    cancerWriter.write(String.format(",%.1f", tpmTotal != null ? tpmTotal : 0));
                }

                cancerWriter.newLine();

                sampleWriter.write(String.format("%s,%s", geneId, geneName));

                final Map<String,Double> sampleTpms = mGeneSampleExpressionData.get(geneId);

                for(final String sampleId : mSampleIds)
                {
                    Double tpm = sampleTpms.get(sampleId);
                    sampleWriter.write(String.format(",%.1f", tpm != null ? tpm : 0));
                }

                sampleWriter.newLine();
            }

            closeBufferedWriter(cancerWriter);
            closeBufferedWriter(sampleWriter);
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write ref RNA gene expression output: {}", e.toString());
        }
    }
}
