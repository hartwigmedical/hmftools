package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.FragmentSizeCalcs.FL_LENGTH;
import static com.hartwig.hmftools.svtools.rna_expression.FragmentSizeCalcs.FL_FREQUENCY;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.GC_ALT;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.GC_CHIMERIC;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.GC_DUPLICATES;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.GC_UNSPLICED;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.GC_READ_THROUGH;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.GC_TOTAL;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.GC_TRANS_SUPPORTING;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TRANS_COUNT;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.UNIQUE_TRANS_COUNT;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ResultsWriter
{
    private final RnaExpConfig mConfig;

    private String mSampledId;

    private BufferedWriter mGeneDataWriter;
    private BufferedWriter mTransDataWriter;
    private BufferedWriter mExonDataWriter;

    private static final Logger LOGGER = LogManager.getLogger(RnaExpression.class);

    public ResultsWriter(final RnaExpConfig config)
    {
        mConfig = config;
        mSampledId = "";

        mGeneDataWriter = null;
        mTransDataWriter = null;
        mExonDataWriter = null;
    }

    public void setSampleId(final String sampleId) { mSampledId = sampleId; }

    public void close()
    {
        closeBufferedWriter(mGeneDataWriter);
        closeBufferedWriter(mTransDataWriter);
        closeBufferedWriter(mExonDataWriter);
    }

    public void writeGeneData(final GeneReadData geneReadData)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mGeneDataWriter == null)
            {
                final String outputFileName = mConfig.OutputDir + "RNA_EXP_GENE_DATA.csv";

                mGeneDataWriter = createBufferedWriter(outputFileName, false);
                mGeneDataWriter.write("SampleId,GeneId,GeneName,Chromosome,GeneLength,TransCount");
                mGeneDataWriter.write(",TotalFragments,SupportingTrans,Alt,Intronic,ReadThrough,Chimeric,Duplicates");
                mGeneDataWriter.newLine();
            }

            final EnsemblGeneData geneData = geneReadData.GeneData;

            mGeneDataWriter.write(String.format("%s,%s,%s,%s,%d,%d",
                    mSampledId, geneData.GeneId, geneData.GeneName, geneData.Chromosome,
                    geneData.GeneEnd - geneData.GeneStart, geneReadData.getTranscripts().size()));

            final int[] fragmentCounts = geneReadData.getCounts();

            mGeneDataWriter.write(String.format(",%d,%d,%d,%d,%d,%d,%d",
                    fragmentCounts[GC_TOTAL], fragmentCounts[GC_TRANS_SUPPORTING], fragmentCounts[GC_ALT],
                    fragmentCounts[GC_UNSPLICED], fragmentCounts[GC_READ_THROUGH],
                    fragmentCounts[GC_CHIMERIC], fragmentCounts[GC_DUPLICATES]));

            mGeneDataWriter.newLine();

        }
        catch(IOException e)
        {
            LOGGER.error("failed to write gene data file: {}", e.toString());
        }
    }

    public void writeTranscriptResults(final GeneReadData geneReadData, final TranscriptResults transResults)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mTransDataWriter == null)
            {
                final String outputFileName = mConfig.OutputDir + "RNA_EXP_TRANS_DATA.csv";

                mTransDataWriter = createBufferedWriter(outputFileName, false);
                mTransDataWriter.write("SampleId,GeneId,GeneName,TransId,Canonical,ExonCount");
                mTransDataWriter.write(",ExonsMatched,ExonicBases,ExonicCoverage,SpliceJuncSupported");
                mTransDataWriter.write(",ShortFragments,ShortUniqueFragments,LongFragments,LongUniqueFragments");
                mTransDataWriter.write(",SpliceJuncFragments,UniqueSpliceJuncFragments");
                mTransDataWriter.newLine();
            }

            final TranscriptData transData = transResults.trans();

            mTransDataWriter.write(String.format("%s,%s,%s,%s,%s,%d",
                    mSampledId, geneReadData.GeneData.GeneId, geneReadData.GeneData.GeneName,
                    transData.TransName, transData.IsCanonical, transData.exons().size()));

            mTransDataWriter.write(String.format(",%d,%d,%d,%d",
                    transResults.exonsFound(), transResults.exonicBases(),
                    transResults.exonicBaseCoverage(), transResults.spliceJunctionsSupported()));

            mTransDataWriter.write(String.format(",%d,%d,%d,%d,%d,%d",
                    transResults.shortSupportingFragments(), transResults.shortUniqueFragments(),
                    transResults.longSupportingFragments(), transResults.longUniqueFragments(),
                    transResults.spliceJunctionFragments(), transResults.spliceJunctionUniqueFragments()));

            mTransDataWriter.newLine();

        }
        catch(IOException e)
        {
            LOGGER.error("failed to write transcripts data file: {}", e.toString());
        }
    }

    public void writeExonData(final GeneReadData geneReadData, final TranscriptData transData)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mExonDataWriter == null)
            {
                final String outputFileName = mConfig.OutputDir + "RNA_EXP_EXON_DATA.csv";

                mExonDataWriter = createBufferedWriter(outputFileName, false);
                mExonDataWriter.write("SampleId,GeneId,GeneName,TransId,ExonRank,ExonStart,ExonEnd");
                mExonDataWriter.write(",TotalCoverage,AvgDepth,Fragments,UniqueFragments");
                mExonDataWriter.write(",SpliceJuncStart,SpliceJuncEnd,UniqueSpliceJuncStart,UniqueSpliceJuncEnd");
                mExonDataWriter.newLine();
            }

            final List<ExonData> exons = transData.exons();

            for(int i = 0; i < exons.size(); ++i)
            {
                ExonData exon = exons.get(i);

                final RegionReadData exonReadData = geneReadData.findExonRegion(exon.ExonStart, exon.ExonEnd);
                if (exonReadData == null)
                    continue;

                mExonDataWriter.write(String.format("%s,%s,%s,%s,%d,%d,%d",
                        mSampledId, geneReadData.GeneData.GeneId, geneReadData.GeneData.GeneName,
                        transData.TransName, exon.ExonRank, exon.ExonStart, exon.ExonEnd));

                int[] matchCounts = exonReadData.getTranscriptReadCount(transData.TransName);
                int[] startSjCounts = exonReadData.getTranscriptJunctionMatchCount(transData.TransName, SE_START);
                int[] endSjCounts = exonReadData.getTranscriptJunctionMatchCount(transData.TransName, SE_END);

                mExonDataWriter.write(String.format(",%d,%.0f,%d,%d,%d,%d,%d,%d",
                        exonReadData.baseCoverage(1), exonReadData.averageDepth(),
                        matchCounts[TRANS_COUNT], matchCounts[UNIQUE_TRANS_COUNT],
                        startSjCounts[TRANS_COUNT], endSjCounts[TRANS_COUNT],
                        startSjCounts[UNIQUE_TRANS_COUNT], endSjCounts[UNIQUE_TRANS_COUNT]));

                /* match types are not recorded per transcript
                mExonDataWriter.write(String.format(",%d,%.0f,%d,%d,%d,%d,%d",
                        exonReadData.matchedReadCount(MATCH_TYPE_WITHIN_EXON),
                        exonReadData.matchedReadCount(MATCH_TYPE_EXON_BOUNDARY) + exonReadData.matchedReadCount(MATCH_TYPE_EXON_MATCH),
                        exonReadData.matchedReadCount(MATCH_TYPE_UNSPLICED));
                 */

                mExonDataWriter.newLine();
            }
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write exon expression file: {}", e.toString());
        }
    }

}
