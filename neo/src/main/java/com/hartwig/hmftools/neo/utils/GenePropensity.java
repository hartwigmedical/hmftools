package com.hartwig.hmftools.neo.utils;

import static com.hartwig.hmftools.common.gene.TranscriptUtils.codingBaseLength;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.formFilename;
import static com.hartwig.hmftools.neo.bind.BindConstants.STRONG_BINDER_LIKELIHOOD;
import static com.hartwig.hmftools.neo.utils.PeptideExpressionData.SOURCE_VALIDATION;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class GenePropensity
{
    private final ExpressionDistribution mExpressionDistribution;
    private final EnsemblDataCache mEnsemblDataCache;

    private final Map<String,TranscriptData> mTransDataMap; // transName to geneId
    private final Map<String, TransTpmLikelihood> mTransLikelihoods; // keyed by transName

    public GenePropensity(final ExpressionDistribution expressionDistribution, final String ensemblDir)
    {
        mExpressionDistribution = expressionDistribution;

        mEnsemblDataCache = new EnsemblDataCache(ensemblDir, RefGenomeVersion.V37);

        mEnsemblDataCache.setRequiredData(true, false, false, false);
        mEnsemblDataCache.load(false);

        mTransDataMap = Maps.newHashMap();
        mTransLikelihoods = Maps.newHashMap();

        for(Map.Entry<String,List<TranscriptData>> entry : mEnsemblDataCache.getTranscriptDataMap().entrySet())
        {
            for(TranscriptData transData : entry.getValue())
            {
                mTransDataMap.put(transData.TransName, transData);
            }
        }
    }

    public void process(final PeptideExpressionData pepExpData)
    {
        if(!pepExpData.hasTpm())
            return;

        double peptideTpm = pepExpData.tpm();

        boolean validationPeptide = pepExpData.Source.equals(SOURCE_VALIDATION);

        if(validationPeptide && pepExpData.LikelihoodRank > STRONG_BINDER_LIKELIHOOD)
            return;

        for(String transName : pepExpData.Transcripts)
        {
            final TranscriptData transData = mTransDataMap.get(transName);

            if(transData == null)
                continue;

            Double transTpm = pepExpData.getTranscriptTpms().get(transName);

            if(transTpm == null)
                continue;

            double expressionLikelihood = mExpressionDistribution.calcLikelihood(transTpm);

            if(expressionLikelihood < 0 || expressionLikelihood > 1)
            {
                NE_LOGGER.error(String.format("peptide(%s) trans(%s) tpm(%.4f) has invalid expression rate(%.4f)",
                        pepExpData, transName, transTpm, expressionLikelihood));
                return;
            }

            TransTpmLikelihood transTpmLikelihood = mTransLikelihoods.get(transName);

            if(transTpmLikelihood == null)
            {
                transTpmLikelihood = new TransTpmLikelihood(transData);
                transTpmLikelihood.ExonCount = transData.exons().size();
                transTpmLikelihood.CodingBases = codingBaseLength(transData);

                mTransLikelihoods.put(transName, transTpmLikelihood);
            }

            if(peptideTpm <= 0)
                continue;

            double transTpmFraction = transTpm / peptideTpm;

            if(validationPeptide)
            {
                ++transTpmLikelihood.ValidationCount;
                transTpmLikelihood.ValidationApp += transTpmFraction;
            }
            else
            {
                ++transTpmLikelihood.StrongBindersCount;
                transTpmLikelihood.StrongBindersApp += transTpmFraction;
                transTpmLikelihood.StrongBindersTpmLikelihoodTotalRaw += expressionLikelihood;
                transTpmLikelihood.StrongBindersTpmLikelihoodTotalApp += expressionLikelihood * transTpmFraction;
            }
        }
    }

    public void writeResults(final String outputDir, final String outputId)
    {
        final String filename = formFilename(outputDir, "trans_propensity", outputId);

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("GeneId,GeneName,TransName,ValidationCount,ValidationnApp,StrongBindersCount,StrongBindersApp");
            writer.write(",SbTpmLikelihoodTotalRaw,SbTpmLikelihoodTotalApp,ExonCount,CodingBases");
            writer.newLine();

            for(TransTpmLikelihood transTpmLikelihood : mTransLikelihoods.values())
            {
                String geneName = mEnsemblDataCache.getGeneDataById(transTpmLikelihood.TransData.GeneId).GeneName;
                writer.write(String.format("%s,%s,%s,%d,%.2f,%d,%.2f,%.4f,%.4f,%d,%d",
                        transTpmLikelihood.TransData.GeneId, geneName, transTpmLikelihood.TransData.TransName,
                        transTpmLikelihood.ValidationCount, transTpmLikelihood.ValidationApp,
                        transTpmLikelihood.StrongBindersCount, transTpmLikelihood.StrongBindersApp,
                        transTpmLikelihood.StrongBindersTpmLikelihoodTotalRaw, transTpmLikelihood.StrongBindersTpmLikelihoodTotalApp,
                        transTpmLikelihood.ExonCount, transTpmLikelihood.CodingBases));

                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write gene propensity file: {}", e.toString());
        }
    }

    private class TransTpmLikelihood
    {
        public final TranscriptData TransData;
        public int ValidationCount;
        public double ValidationApp;
        public int StrongBindersCount;
        public double StrongBindersApp;
        public double StrongBindersTpmLikelihoodTotalRaw;
        public double StrongBindersTpmLikelihoodTotalApp;

        public int CodingBases;
        public int ExonCount;

        public TransTpmLikelihood(final TranscriptData transData)
        {
            TransData = transData;
            ValidationCount = 0;
            ValidationApp = 0;
            StrongBindersCount = 0;
            StrongBindersApp = 0;
            StrongBindersTpmLikelihoodTotalRaw = 0;
            StrongBindersTpmLikelihoodTotalApp = 0;
            CodingBases = 0;
            ExonCount = 0;
        }
    }
}
