package com.hartwig.hmftools.neo.scorer;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.NeoCommon.transcriptsToStr;
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
        mNeoWriter = initNeoepitopeWriter();
    }

    public void close()
    {
        closeBufferedWriter(mNeoWriter);
        closeBufferedWriter(mPeptideWriter);
    }

    private BufferedWriter initPeptideWriter()
    {
        if(!mConfig.WriteTypes.contains(OutputType.ALLELE_PEPTIDE))
            return null;

        try
        {
            final String outputFileName = mConfig.formFilename("peptide_scores");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            if(mConfig.Samples.size() > 1)
                writer.write("SampleId,");

            writer.write("NeIds,VariantType,VariantInfo,Gene,Allele,Peptide,FlankUp,FlankDown");
            writer.write(",Score,Rank,Likelihood,LikelihoodRank,ExpLikelihood,ExpLikelihoodRank,EffectiveTpm");
            writer.write(",RecogSim,OtherAlleleRecogSim,AllelCN,AlleleDisrupted");
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
                if(peptideScoreData.written())
                    continue;

                peptideScoreData.setWritten();

                for(BindData bindData : peptideScoreData.alleleScoreData())
                {
                    if(!logPeptide(bindData))
                        continue;

                    AlleleCoverage alleleCoverage = alleleCoverages.stream().filter(x -> x.Allele.equals(bindData.Allele)).findFirst().orElse(null);

                    if(mConfig.Samples.size() > 1)
                    {
                        mPeptideWriter.write(format("%s,",sampleId));
                    }

                    mPeptideWriter.write(format("%s,%s,%s,%s,%s,%s,%s,%s",
                            peptideScoreData.neoIdsStr(), neoData.VariantType, neoData.VariantInfo, neoData.GeneName,
                            bindData.Allele, bindData.Peptide, bindData.UpFlank, bindData.DownFlank));

                    mPeptideWriter.write(format(",%.4f,%.6f,%.6f,%.6f,%.6f,%.6f,%4.3e,%.1f,%.1f",
                            bindData.score(), bindData.rankPercentile(), bindData.likelihood(), bindData.likelihoodRank(),
                            bindData.expressionLikelihood(), bindData.expressionLikelihoodRank(), peptideScoreData.effectiveTpm(),
                            bindData.recognitionSimilarity(), bindData.otherAlleleRecognitionSimilarity()));

                    mPeptideWriter.write(format(",%.2f,%s", alleleCoverage.CopyNumber, alleleCoverage.isLost()));

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

    private BufferedWriter initNeoepitopeWriter()
    {
        if(!mConfig.WriteTypes.contains(OutputType.NEOEPITOPE))
            return null;

        try
        {
            final String outputFileName = mConfig.formFilename("neoepitope");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            if(mConfig.Samples.size() > 1)
                writer.write("SampleId,");

            writer.write("NeId,VariantType,VariantInfo,GeneName,UpAminoAcids,NovelAminoAcids,DownAminoAcids,PeptideCount");
            writer.write(",TmpSource,RnaFrags,RnaDepth,TpmUp,TpmDown,ExpectedTpm,RawEffectiveTpm,EffectiveTpm");
            writer.write(",TpmCancerUp,TpmCancerDown,TpmPanCancerUp,TpmPanCancerDown");
            writer.write(",NmdMin,NmdMax,CodingBasesLengthMin,CodingBasesLengthMax,FusedIntronLength,SkippedDonors,SkippedAcceptors");
            writer.write(",TranscriptsUp,TranscriptsDown,VariantCopyNumber,CopyNumber,SubclonalLikelihood");


            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to create neoepitope writer: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeNeoData(final String sampleId, final TpmSource tmpSource, final NeoEpitopeData neoData)
    {
        if(mNeoWriter == null)
            return;

        try
        {
            if(mConfig.Samples.size() > 1)
                mNeoWriter.write(format("%s,", sampleId));

            mNeoWriter.write(format("%d,%s,%s,%s,%s,%s,%s",
                    neoData.Id, neoData.VariantType, neoData.VariantInfo, neoData.GeneName,
                    neoData.UpAminoAcids, neoData.NovelAminoAcids, neoData.DownAminoAcids));

            mNeoWriter.write(format(",%d", neoData.peptides().size()));

            final NeoRnaData rnaData = neoData.RnaData;

            mNeoWriter.write(format(",%s,%d,%.0f,%4.3e,%4.3e,%4.3e,%4.3e,%4.3e",
                    tmpSource, rnaData.fragmentSupport(), rnaData.averageBaseDepth(),
                    rnaData.transExpression()[FS_UP], rnaData.transExpression()[FS_DOWN],
                    neoData.expectedTpm(), neoData.rawEffectiveTpm(), neoData.effectiveTpm()));

            mNeoWriter.write(format(",%4.3e,%4.3e,%4.3e,%4.3e",
                    rnaData.tpmCancer()[FS_UP], rnaData.tpmCancer()[FS_DOWN],
                    rnaData.tpmPanCancer()[FS_UP], rnaData.tpmPanCancer()[FS_DOWN]));

            mNeoWriter.write(format(",%d,%d,%d,%d,%d,%d,%d",
                    neoData.NmdBases[FS_UP], neoData.NmdBases[FS_DOWN], neoData.CodingBasesLength[FS_UP], neoData.CodingBasesLength[FS_DOWN],
                    neoData.FusedIntronLength, neoData.SkippedAcceptorsDonors[FS_UP], neoData.SkippedAcceptorsDonors[FS_DOWN]));

            mNeoWriter.write(format(",%s,%s,%.4f,%.4f,%.4f",
                    transcriptsToStr(neoData.Transcripts[FS_UP]), transcriptsToStr(neoData.Transcripts[FS_DOWN]),
                    neoData.VariantCopyNumber, neoData.CopyNumber, neoData.SubclonalLikelihood));

            mNeoWriter.newLine();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write neoepitope data: {}", e.toString());
        }
    }
}
