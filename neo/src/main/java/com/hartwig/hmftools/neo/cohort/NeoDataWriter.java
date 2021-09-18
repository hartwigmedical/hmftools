package com.hartwig.hmftools.neo.cohort;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindScorer.INVALID_CALC;

import java.io.BufferedWriter;
import java.io.IOException;

import com.hartwig.hmftools.neo.bind.BindData;

public class NeoDataWriter
{
    private final NeoScorerConfig mConfig;

    private BufferedWriter mNeoWriter;
    private BufferedWriter mPeptideWriter;

    public NeoDataWriter(final NeoScorerConfig config)
    {
        mConfig = config;

        mPeptideWriter = initPeptideWriter();
        mNeoWriter = null;
        // mNeoWriter = initNeoepitopeWriter();
    }

    public void close()
    {
        closeBufferedWriter(mNeoWriter);
        closeBufferedWriter(mPeptideWriter);
    }

    private BufferedWriter initPeptideWriter()
    {
        try
        {
            final String outputFileName = mConfig.formFilename("allele_peptide_scores");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("SampleId,NeId,VarType,VarInfo,Gene,Allele,Peptide");
            writer.write(",Score,Rank,Likelihood,LikelihoodRank,ExpLikelihood,ExpLikelihoodRank,RecogSim,OtherAlleleRecogSim");
            writer.write(",AllelCN,AlleleDisrupted");
            writer.write(",TpmUp,TpmDown,TpmCancer,TpmCohort,RnaFrags,RnaDepth");
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to create peptide writer: {}", e.toString());
            return null;
        }
    }

    private boolean logPeptide(final BindData bindData)
    {
        if(bindData.likelihoodRank() != INVALID_CALC && bindData.likelihoodRank() <= mConfig.LikelihoodThreshold)
            return true;

        if(bindData.recognitionSimilarity() != INVALID_CALC && bindData.recognitionSimilarity() <= mConfig.SimilarityThreshold)
            return true;

        if(bindData.otherAlleleRecognitionSimilarity() != INVALID_CALC && bindData.otherAlleleRecognitionSimilarity() <= mConfig.SimilarityThreshold)
            return true;

        return false;
    }

    public synchronized void writePeptideData(
            final String sampleId, final NeoEpitopeData neoData, final NeoPredictionData predData, final AlleleCoverage alleleCoverage)
    {
        if(mPeptideWriter == null)
            return;

        try
        {
            for(BindData bindData : predData.getPeptidePredictions(alleleCoverage.Allele))
            {

                if(!logPeptide(bindData))
                    continue;

                mPeptideWriter.write(String.format("%s,%d,%s,%s,%s,%s,%s",
                        sampleId, neoData.Id, neoData.VariantType, neoData.VariantInfo, neoData.GeneName,
                        bindData.Allele, bindData.Peptide));

                mPeptideWriter.write(String.format(",%.4f,%.6f,%.6f,%.6f,%.6f,%.6f,%.1f,%.1f",
                        bindData.score(), bindData.rankPercentile(), bindData.likelihood(), bindData.likelihoodRank(),
                        bindData.expressionLikelihood(), bindData.expressionLikelihoodRank(),
                        bindData.recognitionSimilarity(), bindData.otherAlleleRecognitionSimilarity()));

                mPeptideWriter.write(String.format(",%.2f,%s", alleleCoverage.CopyNumber, alleleCoverage.isLost()));

                mPeptideWriter.write(String.format(",%4.3e,%4.3e,%4.3e,%4.3e,%d,%.0f",
                        neoData.TransExpression[FS_UP], neoData.TransExpression[FS_DOWN], neoData.TpmCancer, neoData.TpmCohort,
                        neoData.RnaNovelFragments, (neoData.RnaBaseDepth[SE_START] + neoData.RnaBaseDepth[SE_END]) * 0.5));

                mPeptideWriter.newLine();
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write peptide data: {}", e.toString());
        }
    }

    /*
    private BufferedWriter initNeoepitopeWriter()
    {
        try
        {
            final String outputFileName = mConfig.formFilename("neoepitope");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("SampleId,NeId,VariantType,VariantInfo,GeneName,AminoAcids");
            writer.write(",Allele,PeptideCount,MaxLikelihood,SumLikelihood");
            writer.write(",AllelCN,AlleleDisrupted");
            writer.write(",TpmSampleUp,TpmSampleDown,TpmCancer,TpmCohort,RnaFrags,RnaDepth");
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to create neoepitope writer: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeNeoData(
            final String sampleId, final NeoEpitopeData neoData, final NeoPredictionData allelePredData, final AlleleCoverage alleleCoverage)
    {
        if(mNeoWriter == null)
            return;

        try
        {
            mNeoWriter.write(String.format("%s,%d,%s,%s,%s,%s",
                    sampleId, neoData.Id, neoData.VariantType, neoData.VariantInfo, neoData.GeneName, neoData.AminoAcids));

            mNeoWriter.write(String.format(",%s,%d,%.4f,%.4f",
                    alleleCoverage.Allele, allelePredData.Peptides, allelePredData.MaxLikelihood, allelePredData.SumLikelihood));

            mNeoWriter.write(String.format(",%.2f,%s", alleleCoverage.CopyNumber, alleleCoverage.isLost()));

            mNeoWriter.write(String.format(",%4.3e,%4.3e,%4.3e,%4.3e,%d,%.0f",
                    neoData.TransExpression[FS_UP], neoData.TransExpression[FS_DOWN], neoData.TpmCancer, neoData.TpmCohort,
                    neoData.RnaNovelFragments, (neoData.RnaBaseDepth[SE_START] + neoData.RnaBaseDepth[SE_END]) * 0.5));

            mNeoWriter.newLine();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write neo-epitope data: {}", e.toString());
        }
    }
    */
}
