package com.hartwig.hmftools.neo.scorer;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindScorer.INVALID_CALC;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

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
            writer.write(",AllelCN,AlleleDisrupted,RnaFrags,RnaDepth");
            writer.write(",TpmUp,TpmDown,TpmDownTotal,ExpectedTpm,EffectiveTpm,RawEffectiveTpm");
            writer.write(",TpmCancerUp,TpmCancerDown,TpmPanCancerUp,TpmPanCancerDown");
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to create peptide writer: {}", e.toString());
            return null;
        }
    }

    public synchronized void writePeptideData(final String sampleId, final NeoEpitopeData neoData, final List<AlleleCoverage> alleleCoverages)
    {
        if(mPeptideWriter == null)
            return;

        try
        {
            for(PeptideScoreData peptideScoreData : neoData.peptides())
            {
                double expectedTpm = peptideScoreData.tpmUpAllocation(neoData.Id);

                for(BindData bindData : peptideScoreData.alleleScoreData())
                {
                    if(!logPeptide(bindData))
                        continue;

                    AlleleCoverage alleleCoverage = alleleCoverages.stream().filter(x -> x.Allele.equals(bindData.Allele)).findFirst().orElse(null);

                    mPeptideWriter.write(String.format("%s,%d,%s,%s,%s,%s,%s",
                            sampleId, neoData.Id, neoData.VariantType, neoData.VariantInfo, neoData.GeneName,
                            bindData.Allele, bindData.Peptide));

                    mPeptideWriter.write(String.format(",%.4f,%.6f,%.6f,%.6f,%.6f,%.6f,%.1f,%.1f",
                            bindData.score(), bindData.rankPercentile(), bindData.likelihood(), bindData.likelihoodRank(),
                            bindData.expressionLikelihood(), bindData.expressionLikelihoodRank(),
                            bindData.recognitionSimilarity(), bindData.otherAlleleRecognitionSimilarity()));

                    mPeptideWriter.write(String.format(",%.2f,%s", alleleCoverage.CopyNumber, alleleCoverage.isLost()));

                    final NeoRnaData rnaData = neoData.RnaData;

                    mPeptideWriter.write(String.format(",%d,%.0f,%4.3e,%4.3e,%4.3e,%4.3e,%4.3e,%4.3e,%4.3e,%4.3e,%4.3e,%4.3e",
                            rnaData.fragmentSupport(), rnaData.averageBaseDepth(),
                            rnaData.transExpression()[FS_UP], rnaData.transExpression()[FS_DOWN], peptideScoreData.tpmDownTotal(),
                            expectedTpm, peptideScoreData.effectiveTpm(), peptideScoreData.rawEffectiveTpm(),
                            rnaData.tpmCancer()[FS_UP], rnaData.tpmCancer()[FS_DOWN],
                            rnaData.tpmPanCancer()[FS_UP], rnaData.tpmPanCancer()[FS_DOWN]));

                    mPeptideWriter.newLine();
                }
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write peptide data: {}", e.toString());
        }
    }

    private boolean logPeptide(final BindData bindData)
    {
        if(passesThreshold(bindData.likelihoodRank(), mConfig.LikelihoodThreshold))
            return true;

        if(passesThreshold(bindData.recognitionSimilarity(), mConfig.SimilarityThreshold))
            return true;

        if(passesThreshold(bindData.otherAlleleRecognitionSimilarity(), mConfig.SimilarityThreshold))
            return true;

        return false;
    }

    public static boolean passesThreshold(double value, double threshold)
    {
        return value != INVALID_CALC && (threshold == 0 || threshold > 0 && value <= threshold);
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

        sj.add(String.format("%6.3e", neo.CancerTpmTotal[FS_UP]));
        sj.add(String.format("%6.3e", neo.CohortTpmTotal[FS_UP]));
        sj.add(String.format("%6.3e", neo.CancerTpmTotal[FS_DOWN]));
        sj.add(String.format("%6.3e", neo.CohortTpmTotal[FS_DOWN]));

            mNeoWriter.newLine();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write neo-epitope data: {}", e.toString());
        }
    }
    */
}
