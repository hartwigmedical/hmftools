package com.hartwig.hmftools.neo.cohort;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.cohort.StatusResults.STATUS_MAX;

import java.io.BufferedWriter;
import java.io.IOException;

public class CohortWriters
{
    private final NeoCohortConfig mConfig;

    private BufferedWriter mNeoWriter;
    private BufferedWriter mPeptideWriter;
    private BufferedWriter mSampleWriter;

    public CohortWriters(final NeoCohortConfig config)
    {
        mConfig = config;

        mNeoWriter = null;
        mPeptideWriter = initPeptideWriter();
        mNeoWriter = initNeoepitopeWriter();

        mSampleWriter = null;
        // initialiseSampleWriter();
    }

    public void close()
    {
        closeBufferedWriter(mNeoWriter);
        closeBufferedWriter(mPeptideWriter);

    }

    private BufferedWriter initNeoepitopeWriter()
    {
        try
        {
            final String outputFileName = mConfig.formFilename("neoepitope");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("SampleId,NeId,VariantType,VariantInfo,AminoAcids");
            writer.write(",Allele,PeptideCount,MaxLikelihood,SumLikelihood");
            writer.write(",MinAffinity,SumAffinity,MinPresentationPerc,SumPresentation");
            writer.write(",AllelCN,AlleleDisrupted");
            writer.write(",TpmSample,TpmCancer,TpmCohort,RnaFrags,RnaDepth");
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
            mNeoWriter.write(String.format("%s,%d,%s,%s,%s",
                    sampleId, neoData.Id, neoData.VariantType, neoData.VariantInfo, neoData.AminoAcids));

            mNeoWriter.write(String.format(",%s,%d,%.4f,%.4f,%.1f,%.4f,%.6f,%.4f",
                    alleleCoverage.Allele, allelePredData.Peptides, allelePredData.MaxLikelihood, allelePredData.SumLikelihood,
                    allelePredData.MinAffinity, allelePredData.SumAffinity, allelePredData.MinPresentationPerc, allelePredData.SumPresentation));

            mNeoWriter.write(String.format(",%.2f,%s", alleleCoverage.CopyNumber, alleleCoverage.isLost()));

            mNeoWriter.write(String.format(",%4.3e,%4.3e,%4.3e,%d,%.0f",
                    neoData.GeneExpression, neoData.TpmCancer, neoData.TpmCohort,
                    neoData.RnaNovelFragments, (neoData.RnaBaseDepth[SE_START] + neoData.RnaBaseDepth[SE_END]) * 0.5));

            mNeoWriter.newLine();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write neo-epitope data: {}", e.toString());
        }
    }

    private BufferedWriter initPeptideWriter()
    {
        try
        {
            final String outputFileName = mConfig.formFilename("allele_peptide");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("SampleId,NeId,Allele,Peptide");
            writer.write(",Score,Rank,Likelihood");
            writer.write(",Affinity,AffinityPerc,PresScore,PresPerc");
            writer.write(",AllelCN,AlleleDisrupted");
            writer.write(",TpmCancer,TpmCohort,RnaFrags,RnaDepth");
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to create peptide writer: {}", e.toString());
            return null;
        }
    }

    public synchronized void writePeptideData(
            final String sampleId, final NeoEpitopeData neoData, final BindingPredictionData predData, final AlleleCoverage alleleCoverage)
    {
        if(mPeptideWriter == null)
            return;

        if(predData.likelihood() < mConfig.LikelihoodThreshold)
            return;

        try
        {
            mPeptideWriter.write(String.format("%s,%d,%s,%s", sampleId, neoData.Id, predData.Allele, predData.Peptide));

            mPeptideWriter.write(String.format(",%.2f,%.6f,%.6f", predData.score(), predData.rankPercentile(), predData.likelihood()));

            mPeptideWriter.write(String.format(",%.2f,%.6f,%.4f,%.6f",
                    predData.affinity(), predData.affinityPerc(), predData.presentation(), predData.presentationPerc()));

            mPeptideWriter.write(String.format(",%.2f,%s", alleleCoverage.CopyNumber, alleleCoverage.isLost()));

            mPeptideWriter.write(String.format(",%4.3e,%4.3e,%d,%.0f",
                    neoData.TpmCancer, neoData.TpmCohort,
                    neoData.RnaNovelFragments, (neoData.RnaBaseDepth[SE_START] + neoData.RnaBaseDepth[SE_END]) * 0.5));

            mPeptideWriter.newLine();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write peptide data: {}", e.toString());
        }
    }

    private void initialiseSampleWriter()
    {
        try
        {
            final String outputFileName = mConfig.OutputDir + "NEO_SAMPLE_SUMMARY.csv";

            mNeoWriter = createBufferedWriter(outputFileName, false);
            mNeoWriter.write("SampleId,PeptideCount," + SampleSummary.header());
            mNeoWriter.newLine();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to create sample summary writer: {}", e.toString());
        }
    }

    private void writeSampleSummary(final String sampleId, final SampleSummary sampleSummary)
    {
        try
        {
            mNeoWriter.write(String.format("%s,%d",
                    sampleId, sampleSummary.PeptideCount));

            for(int i = 0; i < STATUS_MAX; ++i)
            {
                final StatusResults results = sampleSummary.Results[i];

                mNeoWriter.write(String.format(",%.4g,%d,%d,%4g,%d",
                        results.AffinityTotal, results.AffinityLowCount, results.AffinityMediumCount,
                        results.PresentationTotal, results.PresentationCount));
            }

            mNeoWriter.newLine();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write neo-epitope data: {}", e.toString());
        }
    }


}
