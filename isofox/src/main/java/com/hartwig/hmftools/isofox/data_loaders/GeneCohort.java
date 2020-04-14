package com.hartwig.hmftools.isofox.data_loaders;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RnaUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.data_loaders.DataLoaderConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.common.RnaUtils.calcPercentileValues;
import static com.hartwig.hmftools.isofox.data_loaders.TransExpressionCohort.convertDistribution;
import static com.hartwig.hmftools.isofox.data_loaders.TransExpressionCohort.roundTPM;
import static com.hartwig.hmftools.isofox.results.GeneResult.FLD_SUPPORTING_TRANS;
import static com.hartwig.hmftools.isofox.results.GeneResult.FLD_UNSPLICED;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_ID;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class GeneCohort
{
    private final DataLoaderConfig mConfig;

    private final Map<String,GeneCohortData> mGeneCohortDataMap;

    private BufferedWriter mGeneDistributionWriter;

    private static final int DISTRIBUTION_SIZE = 101;

    public GeneCohort(final DataLoaderConfig config)
    {
        mConfig = config;
        mGeneCohortDataMap = Maps.newHashMap();

        mGeneDistributionWriter = null;
    }

    public void processGenes()
    {
        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, DataLoadType.GENE, filenames))
            return;

        initialiseWriter();

        // load each sample's alt SJs and consolidate into a single list
        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path genesFile = filenames.get(i);

            loadFile(sampleId, genesFile);
            ISF_LOGGER.debug("{}: sample({}) loaded genes", i, sampleId);
        }

        ISF_LOGGER.info("loaded {} samples gene files", mConfig.SampleData.SampleIds.size());

        writeGeneFpmPercentiles();

        closeBufferedWriter(mGeneDistributionWriter);
    }

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("gene_distribution.csv");
            mGeneDistributionWriter = createBufferedWriter(outputFileName, false);

            mGeneDistributionWriter.write("GeneId,GeneName");

            double distributionSize = DISTRIBUTION_SIZE;

            for(int i = 0; i <= DISTRIBUTION_SIZE; ++i)
            {
                mGeneDistributionWriter.write(String.format(",Pct_%.2f", i/distributionSize));
            }

            mGeneDistributionWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write gene data file: {}", e.toString());
        }
    }

    private void writeGeneFpmPercentiles()
    {
        try
        {
            int sampleCount = mConfig.SampleData.SampleIds.size();

            for(final GeneCohortData geneData : mGeneCohortDataMap.values())
            {
                final double[] percentileValues = new double[DISTRIBUTION_SIZE + 1];

                calcPercentileValues(convertDistribution(geneData.FragsPerMillionValues, sampleCount), percentileValues);

                mGeneDistributionWriter.write(String.format("%s,%s",geneData.GeneId, geneData.GeneName));

                for (int i = 0; i <= DISTRIBUTION_SIZE; ++i)
                {
                    mGeneDistributionWriter.write(String.format(",%6.3e", percentileValues[i]));
                }

                mGeneDistributionWriter.newLine();
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write gene data file: {}", e.toString());
        }
    }

    private void loadFile(final String sampleId, final Path filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(filename);

            final Map<String,Integer> fieldsMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            int geneIdIndex = fieldsMap.get(FLD_GENE_ID);
            int supportTransIndex = fieldsMap.get(FLD_SUPPORTING_TRANS);
            int unsplicedIndex = fieldsMap.get(FLD_UNSPLICED);

            final Map<String,double[]> geneFpmData = Maps.newHashMap();

            for(final String data : lines)
            {
                final String[] items = data.split(DELIMITER);

                final String geneId = items[geneIdIndex];

                boolean excludeGene = !mConfig.RestrictedGeneIds.isEmpty() && !mConfig.RestrictedGeneIds.contains(geneId);

                if(!excludeGene)
                {
                    GeneCohortData geneData = mGeneCohortDataMap.get(geneId);

                    if (geneData == null)
                    {
                        geneData = GeneCohortData.fromCsv(items, fieldsMap);
                        mGeneCohortDataMap.put(geneId, geneData);
                    }
                }

                double[] fpmData = new double[FPM_FPM + 1];
                fpmData[FPM_SUPPORTING] = Integer.parseInt(items[supportTransIndex]);
                fpmData[FPM_UNSPLICED] = Integer.parseInt(items[unsplicedIndex]);

                geneFpmData.put(geneId, fpmData);
            }

            calcAndAddFpmData(sampleId, geneFpmData);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load transcript data file({}): {}", filename.toString(), e.toString());
            return;
        }
    }

    private static final int FPM_SUPPORTING = 0;
    private static final int FPM_UNSPLICED = 1;
    private static final int FPM_FPM = 2;

    private void calcAndAddFpmData(final String sampleId, final Map<String,double[]> geneFpmData)
    {
        long totalFragments = geneFpmData.values().stream().mapToLong(x -> (long)x[FPM_SUPPORTING] + (long)x[FPM_UNSPLICED]).sum();
        double fpmFactor = 1000000.0 / totalFragments;

        boolean roundValues = mConfig.SampleData.SampleIds.size() >= 100;

        for(Map.Entry<String,double[]> entry : geneFpmData.entrySet())
        {
            final String geneId = entry.getKey();
            final double[] fpmData = entry.getValue();
            fpmData[FPM_FPM] = (fpmData[FPM_SUPPORTING] + fpmData[FPM_UNSPLICED]) * fpmFactor;

            GeneCohortData geneData = mGeneCohortDataMap.get(geneId);

            if(geneData == null) // filtered out earlier
                continue;

            geneData.addSampleData(
                    sampleId,
                    roundValues ? roundTPM(fpmData[FPM_FPM], mConfig.TpmRounding) : fpmData[FPM_FPM]);
        }
    }


}
