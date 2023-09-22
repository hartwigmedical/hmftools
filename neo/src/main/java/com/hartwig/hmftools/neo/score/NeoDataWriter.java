package com.hartwig.hmftools.neo.score;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.NeoCommon.transcriptsToStr;
import static com.hartwig.hmftools.neo.bind.BindScorer.INVALID_CALC;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

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

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            if(mConfig.Samples.size() > 1)
                sj.add("SampleId");

            sj.add("NeIds").add("VariantType").add("VariantInfo").add("Gene").add("Allele").add("Peptide").add("FlankUp").add("FlankDown");
            sj.add("Score").add("Rank").add("Likelihood").add("LikelihoodRank").add("ExpLikelihood").add("ExpLikelihoodRank").add("EffectiveTpm");
            sj.add("RecogSim").add("OtherAlleleRecogSim").add("AlleleCN").add("AlleleDisrupted");

            writer.write(sj.toString());
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

                    StringJoiner sj = new StringJoiner(TSV_DELIM);

                    if(mConfig.Samples.size() > 1)
                        sj.add(sampleId);

                    sj.add(peptideScoreData.neoIdsStr());
                    sj.add(neoData.VariantType.toString());
                    sj.add(neoData.VariantInfo);
                    sj.add(neoData.GeneName);
                    sj.add(bindData.Allele);
                    sj.add(bindData.Peptide);
                    sj.add(bindData.UpFlank);
                    sj.add(bindData.DownFlank);

                    sj.add(format("%.4f", bindData.score()));
                    sj.add(format("%.6f", bindData.rankPercentile()));
                    sj.add(format("%.6f", bindData.likelihood()));
                    sj.add(format("%.6f", bindData.likelihoodRank()));
                    sj.add(format("%.6f", bindData.expressionLikelihood()));
                    sj.add(format("%.6f", bindData.expressionLikelihoodRank()));
                    sj.add(format("%4.3e", peptideScoreData.effectiveTpm()));
                    sj.add(format("%.1f", bindData.recognitionSimilarity()));
                    sj.add(format("%.1f", bindData.otherAlleleRecognitionSimilarity()));

                    sj.add(format("%.1f", alleleCoverage.CopyNumber));
                    sj.add(format("%s", alleleCoverage.isLost()));

                    mPeptideWriter.write(sj.toString());
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

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            if(mConfig.Samples.size() > 1)
                sj.add("SampleId");

            sj.add("NeId").add("VariantType").add("VariantInfo").add("GeneName")
                    .add("UpAminoAcids").add("NovelAminoAcids").add("DownAminoAcids")
                    .add("PeptideCount").add("TpmSource").add("RnaFrags").add("RnaDepth").add("TpmUp").add("TpmDown").add("ExpectedTpm")
                    .add("RawEffectiveTpm").add("EffectiveTpm").add("TpmCancerUp").add("TpmCancerDown").add("TpmPanCancerUp").add("TpmPanCancerDown")
                    .add("NmdMin").add("NmdMax").add("CodingBasesLengthMin").add("CodingBasesLengthMax").add("FusedIntronLength")
                    .add("SkippedDonors").add("SkippedAcceptors")
                    .add("TranscriptsUp").add("TranscriptsDown").add("VariantCopyNumber").add("CopyNumber").add("SubclonalLikelihood");

            writer.write(sj.toString());
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
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            if(mConfig.Samples.size() > 1)
                sj.add(sampleId);

            sj.add(String.valueOf(neoData.Id));
            sj.add(neoData.VariantType.toString());
            sj.add(neoData.VariantInfo);
            sj.add(neoData.GeneName);
            sj.add(neoData.UpAminoAcids);
            sj.add(neoData.NovelAminoAcids);
            sj.add(neoData.DownAminoAcids);

            sj.add(String.valueOf(neoData.peptides().size()));

            final NeoRnaData rnaData = neoData.RnaData;

            sj.add(tmpSource.toString());
            sj.add(String.valueOf(rnaData.fragmentSupport()));
            sj.add(format("%.0f", rnaData.averageBaseDepth()));
            sj.add(format("%4.3e", rnaData.transExpression()[FS_UP]));
            sj.add(format("%4.3e", rnaData.transExpression()[FS_DOWN]));
            sj.add(format("%4.3e", neoData.expectedTpm()));
            sj.add(format("%4.3e", neoData.rawEffectiveTpm()));
            sj.add(format("%4.3e", neoData.effectiveTpm()));

            sj.add(format("%4.3e", rnaData.tpmCancer()[FS_UP]));
            sj.add(format("%4.3e", rnaData.tpmCancer()[FS_DOWN]));
            sj.add(format("%4.3e", rnaData.tpmPanCancer()[FS_UP]));
            sj.add(format("%4.3e", rnaData.tpmPanCancer()[FS_DOWN]));

            sj.add(String.valueOf(neoData.NmdBases[FS_UP]));
            sj.add(String.valueOf(neoData.NmdBases[FS_DOWN]));
            sj.add(String.valueOf(neoData.CodingBasesLength[FS_UP]));
            sj.add(String.valueOf(neoData.CodingBasesLength[FS_DOWN]));
            sj.add(String.valueOf(neoData.FusedIntronLength));
            sj.add(String.valueOf(neoData.SkippedAcceptorsDonors[FS_UP]));
            sj.add(String.valueOf(neoData.SkippedAcceptorsDonors[FS_DOWN]));

            sj.add(transcriptsToStr(neoData.Transcripts[FS_UP]));
            sj.add(transcriptsToStr(neoData.Transcripts[FS_DOWN]));
            sj.add(format("%.4f", neoData.VariantCopyNumber));
            sj.add(format("%.4f", neoData.CopyNumber));
            sj.add(format("%.4f", neoData.SubclonalLikelihood));

            mNeoWriter.write(sj.toString());
            mNeoWriter.newLine();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write neoepitope data: {}", e.toString());
        }
    }
}
