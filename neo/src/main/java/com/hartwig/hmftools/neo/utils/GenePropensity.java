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
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.compress.utils.Lists;

public class GenePropensity
{
    private final ExpressionDistribution mExpressionDistribution;
    private final EnsemblDataCache mEnsemblDataCache;

    private final Map<String,TranscriptData> mTransDataMap; // transName to geneId
    private final Map<String, GeneTpmLikelihood> mGeneLikelihoods; // keyed by gendId

    public GenePropensity(final ExpressionDistribution expressionDistribution, final String ensemblDir)
    {
        mExpressionDistribution = expressionDistribution;

        mEnsemblDataCache = new EnsemblDataCache(ensemblDir, RefGenomeVersion.V37);

        mEnsemblDataCache.setRequiredData(true, false, false, false);
        mEnsemblDataCache.load(false);

        mTransDataMap = Maps.newHashMap();
        mGeneLikelihoods = Maps.newHashMap();

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

        boolean validationPeptide = pepExpData.Source.equals(SOURCE_VALIDATION);

        if(validationPeptide && pepExpData.LikelihoodRank > STRONG_BINDER_LIKELIHOOD)
            return;

        // organise transcripts into lists by gene
        Map<String,List<TranscriptData>> geneTransMap = Maps.newHashMap();

        for(String transName : pepExpData.Transcripts)
        {
            final TranscriptData transData = mTransDataMap.get(transName);

            if(transData == null)
                continue;

            List<TranscriptData> transList = geneTransMap.get(transData.GeneId);

            if(transList == null)
            {
                transList = Lists.newArrayList();
                geneTransMap.put(transData.GeneId, transList);
            }

            transList.add(transData);
        }

        for(Map.Entry<String,List<TranscriptData>> entry : geneTransMap.entrySet())
        {
            // sum up TPM for this gene then get expression likelihood
            String geneId = entry.getKey();
            List<TranscriptData> transList = entry.getValue();

            double geneTpm = transList.stream().map(x -> pepExpData.getTranscriptTpms().get(x.TransName))
                    .filter(x -> x != null).mapToDouble(x -> x).sum();

            double expressionLikelihood = mExpressionDistribution.calcLikelihood(geneTpm);

            if(expressionLikelihood < 0 || expressionLikelihood > 1)
            {
                NE_LOGGER.error(String.format("peptide(%s) tpm(%.4f) has invalid expression rate(%.4f)",
                        pepExpData, geneTpm, expressionLikelihood));
                return;
            }

            GeneTpmLikelihood geneTpmLikelihood = mGeneLikelihoods.get(geneId);

            if(geneTpmLikelihood == null)
            {
                geneTpmLikelihood = new GeneTpmLikelihood();

                // find canonical
                TranscriptData canonicalTrans = transList.stream().filter(x -> x.IsCanonical).findFirst().orElse(null);

                if(canonicalTrans == null)
                    canonicalTrans = mEnsemblDataCache.getTranscriptData(geneId, "");

                if(canonicalTrans != null)
                {
                    geneTpmLikelihood.ExonCount = canonicalTrans.exons().size();
                    geneTpmLikelihood.CodingBases = codingBaseLength(canonicalTrans);
                }

                mGeneLikelihoods.put(geneId, geneTpmLikelihood);
            }

            if(validationPeptide)
            {
                ++geneTpmLikelihood.ValidationCount;
            }
            else
            {
                ++geneTpmLikelihood.StrongBindersCount;
                geneTpmLikelihood.StrongBindersTpmLikelihoodTotal += expressionLikelihood;
            }
        }
    }

    public void writeResults(final String outputDir, final String outputId)
    {
        final String filename = formFilename(outputDir, "gene_propensity", outputId);

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("GeneId,GeneName,ValidationCount,StrongBindersCount,ExonCount,CodingBases,StrongBindersTpmLikelihoodTotal");
            writer.newLine();

            for(Map.Entry<String,GeneTpmLikelihood> entry : mGeneLikelihoods.entrySet())
            {
                String geneId = entry.getKey();
                GeneTpmLikelihood geneTpmLikelihood = entry.getValue();

                writer.write(String.format("%s,%s,%d,%d,%d,%d,%.4f",
                        geneId, mEnsemblDataCache.getGeneDataById(geneId).GeneName,
                        geneTpmLikelihood.ValidationCount, geneTpmLikelihood.StrongBindersCount,
                        geneTpmLikelihood.ExonCount, geneTpmLikelihood.CodingBases, geneTpmLikelihood.StrongBindersTpmLikelihoodTotal));

                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write gene propensity file: {}", e.toString());
        }
    }

    private class GeneTpmLikelihood
    {
        public int ValidationCount;
        public int StrongBindersCount;
        public double StrongBindersTpmLikelihoodTotal;

        public int CodingBases;
        public int ExonCount;

        public GeneTpmLikelihood()
        {
            ValidationCount = 0;
            StrongBindersCount = 0;
            StrongBindersTpmLikelihoodTotal = 0;
            CodingBases = 0;
            ExonCount = 0;
        }
    }
}
