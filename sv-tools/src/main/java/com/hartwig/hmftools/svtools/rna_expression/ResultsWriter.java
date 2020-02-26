package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.ExpectedTransRates.UNSPLICED_ID;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.ALT;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.CHIMERIC;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.DUPLICATE;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.READ_THROUGH;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.TOTAL;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.TRANS_SUPPORTING;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.UNSPLICED;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.typeAsInt;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TRANS_COUNT;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.UNIQUE_TRANS_COUNT;
import static com.hartwig.hmftools.svtools.rna_expression.TranscriptModel.calcEffectiveLength;

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

    private BufferedWriter mGeneDataWriter;
    private BufferedWriter mTransDataWriter;
    private BufferedWriter mExonDataWriter;
    private BufferedWriter mTransComboWriter;

    // controlled by other components but instantiated once for output synchronosation
    private BufferedWriter mExpRateWriter;
    private BufferedWriter mReadDataWriter;
    private BufferedWriter mAltSpliceJunctionWriter;

    private static final Logger LOGGER = LogManager.getLogger(RnaExpression.class);

    public ResultsWriter(final RnaExpConfig config)
    {
        mConfig = config;

        mGeneDataWriter = null;
        mTransDataWriter = null;
        mExonDataWriter = null;
        mTransComboWriter = null;
        mExpRateWriter = null;
        mReadDataWriter = null;
        mAltSpliceJunctionWriter = null;

        initialiseExternalWriters();
    }

    public void close()
    {
        closeBufferedWriter(mGeneDataWriter);
        closeBufferedWriter(mTransDataWriter);
        closeBufferedWriter(mExonDataWriter);
        closeBufferedWriter(mTransComboWriter);
        closeBufferedWriter(mExpRateWriter);
        closeBufferedWriter(mReadDataWriter);
        closeBufferedWriter(mAltSpliceJunctionWriter);
    }

    private void initialiseExternalWriters()
    {
        if(mConfig.GenerateExpectedRates)
        {
            mExpRateWriter = ExpectedRatesGenerator.createWriter(mConfig);
        }
        else
        {
            if (mConfig.WriteReadData)
                mReadDataWriter = GeneBamReader.createReadDataWriter(mConfig);

            mAltSpliceJunctionWriter = GeneBamReader.createAltSpliceJunctionWriter(mConfig);
        }
    }

    public BufferedWriter getExpRatesWriter() { return mExpRateWriter;}
    public BufferedWriter getAltSpliceJunctionWriter() { return mAltSpliceJunctionWriter;}
    public BufferedWriter getReadDataWriter() { return mReadDataWriter; }

    public synchronized void writeGeneData(final GeneReadData geneReadData)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mGeneDataWriter == null)
            {
                final String outputFileName = mConfig.formOutputFile("gene_data.csv");

                mGeneDataWriter = createBufferedWriter(outputFileName, false);
                mGeneDataWriter.write("GeneId,GeneName,Chromosome,GeneLength,IntronicLength,TransCount");
                mGeneDataWriter.write(",TotalFragments,SupportingTrans,Alt,Unspliced,ReadThrough,Chimeric,Duplicates,UnsplicedAlloc,FitResiduals");
                mGeneDataWriter.newLine();
            }

            final EnsemblGeneData geneData = geneReadData.GeneData;

            long geneLength = geneData.GeneEnd - geneData.GeneStart;

            mGeneDataWriter.write(String.format("%s,%s,%s,%d,%d,%d",
                    geneData.GeneId, geneData.GeneName, geneData.Chromosome, geneLength,
                    geneLength - geneReadData.calcExonicRegionLength(), geneReadData.getTranscripts().size()));

            final int[] fragmentCounts = geneReadData.getCounts();

            mGeneDataWriter.write(String.format(",%d,%d,%d,%d,%d,%d,%d",
                    fragmentCounts[typeAsInt(TOTAL)], fragmentCounts[typeAsInt(TRANS_SUPPORTING)], fragmentCounts[typeAsInt(ALT)],
                    fragmentCounts[typeAsInt(UNSPLICED)], fragmentCounts[typeAsInt(READ_THROUGH)],
                    fragmentCounts[typeAsInt(CHIMERIC)], fragmentCounts[typeAsInt(DUPLICATE)]));

            Double unsplicedAlloc = geneReadData.getTranscriptAllocations().get(UNSPLICED_ID);
            mGeneDataWriter.write(String.format(",%.1f,%.1f",
                    unsplicedAlloc != null && !Double.isNaN(unsplicedAlloc) ? unsplicedAlloc : 0.0,
                    geneReadData.getFitResiduals()));

            mGeneDataWriter.newLine();

        }
        catch(IOException e)
        {
            LOGGER.error("failed to write gene data file: {}", e.toString());
        }
    }

    public synchronized void writeTranscriptResults(final GeneReadData geneReadData, final TranscriptResults transResults)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mTransDataWriter == null)
            {
                final String outputFileName = mConfig.formOutputFile("transcript_data.csv");

                mTransDataWriter = createBufferedWriter(outputFileName, false);
                mTransDataWriter.write("GeneId,GeneName,TransId,TransName,Canonical,ExonCount,EffectiveLength");
                mTransDataWriter.write(",ExonsMatched,ExonicBases,ExonicCoverage,FitAllocation");
                mTransDataWriter.write(",UniqueBases,UniqueBaseCoverage,UniqueBaseAvgDepth");
                mTransDataWriter.write(",SpliceJuncSupported,UniqueSpliceJunc,UniqueSpliceJuncSupported");
                mTransDataWriter.write(",ShortFragments,ShortUniqueFragments,LongFragments,LongUniqueFragments,SpliceJuncFragments,UniqueSpliceJuncFragments");
                mTransDataWriter.newLine();
            }

            final TranscriptData transData = transResults.trans();

            double effectiveLength = calcEffectiveLength(transResults.exonicBases(), mConfig.ExpRateFragmentLengths);

            mTransDataWriter.write(String.format("%s,%s,%d,%s,%s,%d,%.0f",
                    geneReadData.GeneData.GeneId, geneReadData.GeneData.GeneName,
                    transData.TransId, transData.TransName, transData.IsCanonical, transData.exons().size(), effectiveLength));

            double expRateAllocation = geneReadData.getTranscriptAllocation(transData.TransName);

            mTransDataWriter.write(String.format(",%d,%d,%d,%.1f,%d,%d,%.0f",
                    transResults.exonsFound(), transResults.exonicBases(), transResults.exonicBaseCoverage(),
                    expRateAllocation, transResults.uniqueBases(), transResults.uniqueBaseCoverage(), transResults.uniqueBaseAvgDepth()));

            mTransDataWriter.write(String.format(",%d,%d,%d",
                    transResults.spliceJunctionsSupported(), transResults.uniqueSpliceJunctions(), transResults.uniqueSpliceJunctionsSupported()));

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

    public synchronized void writeExonData(final GeneReadData geneReadData, final TranscriptData transData)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mExonDataWriter == null)
            {
                final String outputFileName = mConfig.formOutputFile("exon_data.csv");

                mExonDataWriter = createBufferedWriter(outputFileName, false);
                mExonDataWriter.write("GeneId,GeneName,TransId,TransName,ExonRank,ExonStart,ExonEnd,SharedTrans");
                mExonDataWriter.write(",TotalCoverage,AvgDepth,UniqueBases,UniqueBaseCoverage,UniqueBaseAvgDepth,Fragments,UniqueFragments");
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

                mExonDataWriter.write(String.format("%s,%s,%d,%s",
                        geneReadData.GeneData.GeneId, geneReadData.GeneData.GeneName, transData.TransId, transData.TransName));

                mExonDataWriter.write(String.format(",%d,%d,%d,%d",
                        exon.ExonRank, exon.ExonStart, exon.ExonEnd, exonReadData.getRefRegions().size()));

                int[] matchCounts = exonReadData.getTranscriptReadCount(transData.TransId);
                int[] startSjCounts = exonReadData.getTranscriptJunctionMatchCount(transData.TransId, SE_START);
                int[] endSjCounts = exonReadData.getTranscriptJunctionMatchCount(transData.TransId, SE_END);

                int uniqueBaseTotalDepth = exonReadData.uniqueBaseTotalDepth();
                int uniqueBaseCount = exonReadData.uniqueBaseCount();
                double uniqueAvgDepth = uniqueBaseCount > 0 ? uniqueBaseTotalDepth / (double)uniqueBaseCount : 0;

                mExonDataWriter.write(String.format(",%d,%.0f,%d,%d,%.0f",
                        exonReadData.baseCoverage(1), exonReadData.averageDepth(),
                        uniqueBaseCount, exonReadData.uniqueBaseCoverage(1), uniqueAvgDepth));

                mExonDataWriter.write(String.format(",%d,%d,%d,%d,%d,%d",
                        matchCounts[TRANS_COUNT], matchCounts[UNIQUE_TRANS_COUNT],
                        startSjCounts[TRANS_COUNT], endSjCounts[TRANS_COUNT],
                        startSjCounts[UNIQUE_TRANS_COUNT], endSjCounts[UNIQUE_TRANS_COUNT]));

                mExonDataWriter.newLine();
            }
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write exon expression file: {}", e.toString());
        }
    }

    public synchronized void writeTransComboCounts(
            final GeneReadData geneReadData, final List<String> categories, final double[] counts, final double[] fittedCounts)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mTransComboWriter == null)
            {
                final String outputFileName = mConfig.formOutputFile("category_counts.csv");

                mTransComboWriter = createBufferedWriter(outputFileName, false);
                mTransComboWriter.write("GeneId,GeneName,Category,Count,FitCount");
                mTransComboWriter.newLine();
            }

            for(int i = 0; i < categories.size(); ++i)
            {
                double count = counts[i];
                final String category = categories.get(i);

                mTransComboWriter.write(String.format("%s,%s,%s,%.0f,%.1f",
                        geneReadData.GeneData.GeneId, geneReadData.GeneData.GeneName, category,
                        count, fittedCounts[i]));

                mTransComboWriter.newLine();
            }
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write trans combo data file: {}", e.toString());
        }
    }
}
