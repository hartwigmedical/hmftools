package com.hartwig.hmftools.isofox.expression.cohort;

import static com.hartwig.hmftools.common.sigs.DataUtils.convertList;
import static com.hartwig.hmftools.common.sigs.Percentiles.calcPercentileValues;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.expression.cohort.TransExpressionDistribution.DISTRIBUTION_SIZE;
import static com.hartwig.hmftools.isofox.expression.cohort.TransExpressionDistribution.convertDistribution;
import static com.hartwig.hmftools.isofox.expression.cohort.TransExpressionDistribution.roundTPM;
import static com.hartwig.hmftools.isofox.results.GeneResult.FLD_SPLICED_FRAGS;
import static com.hartwig.hmftools.isofox.results.GeneResult.FLD_UNSPLICED_FRAGS;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_ID;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_TPM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.cohort.CohortAnalysisType;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

public class GeneExpressionDistribution
{
    private final CohortConfig mConfig;

    private final String mGeneRateType;

    private final Map<String,GeneCohortData> mGeneCohortDataMap;

    private BufferedWriter mWriter;

    public static final String GENE_RATE_TPM = "TPM";
    public static final String GENE_RATE_FPM = "FPM";

    public GeneExpressionDistribution(final CohortConfig config)
    {
        mConfig = config;
        mGeneCohortDataMap = Maps.newHashMap();

        mGeneRateType = GENE_RATE_TPM;
        mWriter = null;
    }

    public void processGenes()
    {
        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, CohortAnalysisType.GENE_DISTRIBUTION, filenames))
            return;

        initialiseWriter();

        // load each sample's gene expression records and consolidate into a single list
        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path genesFile = filenames.get(i);

            loadFile(sampleId, genesFile);
            ISF_LOGGER.debug("{}: sample({}) loaded genes", i, sampleId);
        }

        ISF_LOGGER.info("loaded {} samples gene files", mConfig.SampleData.SampleIds.size());

        writeGeneRatePercentiles();

        closeBufferedWriter(mWriter);
    }

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("gene_distribution.csv");
            mWriter = createBufferedWriter(outputFileName, false);

            mWriter.write("GeneId,GeneName");

            for(int i = 0; i < DISTRIBUTION_SIZE; ++i)
            {
                mWriter.write(String.format(",Pct_%d", i));
            }

            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write gene data file: {}", e.toString());
        }
    }

    private void writeGeneRatePercentiles()
    {
        try
        {
            int sampleCount = mConfig.SampleData.SampleIds.size();

            for(final GeneCohortData geneData : mGeneCohortDataMap.values())
            {
                final double[] percentileValues = new double[DISTRIBUTION_SIZE];

                final List<Double> geneRateValues = convertDistribution(geneData.GeneRates, sampleCount);
                final double[] geneRateDataList = convertList(geneRateValues);

                calcPercentileValues(geneRateDataList, percentileValues);

                mWriter.write(String.format("%s,%s",geneData.GeneId, geneData.GeneName));

                for (int i = 0; i < DISTRIBUTION_SIZE; ++i)
                {
                    mWriter.write(String.format(",%6.3e", percentileValues[i]));
                }

                mWriter.newLine();
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
            boolean roundValues = mConfig.SampleData.SampleIds.size() >= 100;

            final List<String> lines = Files.readAllLines(filename);

            final Map<String,Integer> fieldsMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            int geneIdIndex = fieldsMap.get(FLD_GENE_ID);
            int splicedIndex = fieldsMap.get(FLD_SPLICED_FRAGS);
            int unsplicedIndex = fieldsMap.get(FLD_UNSPLICED_FRAGS);
            int tpmIndex = fieldsMap.get(FLD_TPM);

            final Map<String,double[]> geneFpmData = Maps.newHashMap();

            for(final String data : lines)
            {
                final String[] items = data.split(DELIMITER);

                final String geneId = items[geneIdIndex];

                boolean excludeGene = !mConfig.RestrictedGeneIds.isEmpty() && !mConfig.RestrictedGeneIds.contains(geneId);

                GeneCohortData geneData = null;

                if (!excludeGene)
                {
                    geneData = mGeneCohortDataMap.get(geneId);

                    if (geneData == null)
                    {
                        geneData = GeneCohortData.fromCsv(items, fieldsMap);
                        mGeneCohortDataMap.put(geneId, geneData);
                    }
                }

                if (mGeneRateType.equals(GENE_RATE_FPM))
                {
                    double[] fpmData = new double[FPM_FPM + 1];
                    fpmData[FPM_SUPPORTING] = Integer.parseInt(items[splicedIndex]);
                    fpmData[FPM_UNSPLICED] = Integer.parseInt(items[unsplicedIndex]);

                    geneFpmData.put(geneId, fpmData);
                }
                else if (geneData != null)
                {
                    double tpm = Double.parseDouble(items[tpmIndex]);
                    geneData.addSampleData(sampleId, roundValues ? roundTPM(tpm, mConfig.Expression.TpmRounding) : tpm);
                }
            }

            calcAndAddFpmData(sampleId, geneFpmData);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load gene data file({}): {}", filename.toString(), e.toString());
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
                    roundValues ? roundTPM(fpmData[FPM_FPM], mConfig.Expression.TpmRounding) : fpmData[FPM_FPM]);
        }
    }

}
