package com.hartwig.hmftools.neo.missense;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.neo.bind.BindData;

public class MissenseWriter
{
    private final BufferedWriter mWriter;
    private final boolean mRunScoring;
    private final MissenseConfig mConfig;

    public MissenseWriter(final MissenseConfig config, boolean runScoring)
    {
        mRunScoring = runScoring;
        mConfig = config;

        String filePrefix = mRunScoring ? "missense_peptide_scores" : "missense_peptides";
        mWriter = initialiseWriter(formOutputFile(filePrefix));
    }

    private String formOutputFile(final String fileType)
    {
        String filePrefix = mConfig.OutputDir + fileType;

        if(mConfig.OutputId != null)
            filePrefix += "." + mConfig.OutputId;

        return filePrefix + TSV_EXTENSION;
    }

    private BufferedWriter initialiseWriter(final String outputFile)
    {
        NE_LOGGER.info("writing results to {}", outputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("GeneId\tGeneName\tTransName\tPosition\tCodonIndex\tContext\tRef\tAlt\tPeptide\tUpFlank\tDownFlank");

            if(mRunScoring)
                writer.write("\tAllele\tScore\tRank\tLikelihood\tLikelihoodRank\tExpLikelihood\tExpLikelihoodRank");

            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    public synchronized void writePeptideData(final List<MissensePeptide> peptideDataList, final List<BindData> bindDataList)
    {
        try
        {
            for(int i = 0; i < peptideDataList.size(); ++i)
            {
                MissensePeptide peptideData = peptideDataList.get(i);

                mWriter.write(String.format("%s\t%s\t%s",
                        peptideData.GeneId, peptideData.GeneName, peptideData.TransName));

                mWriter.write(String.format("\t%d\t%d\t%s\t%c\t%c",
                        peptideData.Position, peptideData.CodonIndex, peptideData.Context, peptideData.RefBase, peptideData.AltBase));

                mWriter.write(String.format("\t%s\t%s\t%s",
                        peptideData.Peptide, peptideData.UpFlank, peptideData.DownFlank));

                if(mRunScoring)
                {
                    BindData bindData = bindDataList.get(i);

                    mWriter.write(String.format("\t%s\t%.4f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f",
                            bindData.Allele, bindData.score(), bindData.rankPercentile(), bindData.likelihood(), bindData.likelihoodRank(),
                            bindData.expressionLikelihood(), bindData.expressionLikelihoodRank()));
                }

                mWriter.newLine();
            }
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write peptide data: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mWriter);
    }
}
