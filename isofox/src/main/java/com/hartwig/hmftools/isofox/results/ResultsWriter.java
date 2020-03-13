package com.hartwig.hmftools.isofox.results;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.GeneCollection.TRANS_COUNT;
import static com.hartwig.hmftools.isofox.common.GeneCollection.UNIQUE_TRANS_COUNT;
import static com.hartwig.hmftools.isofox.common.RegionReadData.findExonRegion;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.isofox.GeneBamReader;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.FragmentSizeCalcs;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator;
import com.hartwig.hmftools.isofox.gc.GcRatioCounts;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunctionFinder;
import com.hartwig.hmftools.isofox.novel.RetainedIntronFinder;

public class ResultsWriter
{
    private final IsofoxConfig mConfig;

    private BufferedWriter mGeneDataWriter;
    private BufferedWriter mTransDataWriter;
    private BufferedWriter mExonDataWriter;
    private BufferedWriter mTransComboWriter;

    // controlled by other components but instantiated once for output synchronosation
    private BufferedWriter mExpRateWriter;
    private BufferedWriter mReadDataWriter;
    private BufferedWriter mAltSpliceJunctionWriter;
    private BufferedWriter mFragLengthWriter;
    private BufferedWriter mReadGcRatioWriter;
    private BufferedWriter mRetainedIntronWriter;

    public ResultsWriter(final IsofoxConfig config)
    {
        mConfig = config;

        mGeneDataWriter = null;
        mTransDataWriter = null;
        mExonDataWriter = null;
        mTransComboWriter = null;
        mExpRateWriter = null;
        mReadDataWriter = null;
        mAltSpliceJunctionWriter = null;
        mFragLengthWriter = null;
        mReadGcRatioWriter = null;
        mRetainedIntronWriter = null;

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
        closeBufferedWriter(mFragLengthWriter);
        closeBufferedWriter(mReadGcRatioWriter);
        closeBufferedWriter(mRetainedIntronWriter);
    }

    private void initialiseExternalWriters()
    {
        if(mConfig.WriteExpectedRates)
        {
            mExpRateWriter = ExpectedRatesGenerator.createWriter(mConfig);
        }
        else
        {
            if (mConfig.WriteReadData)
                mReadDataWriter = GeneBamReader.createReadDataWriter(mConfig);

            if(mConfig.WriteTransData)
            {
                mAltSpliceJunctionWriter = AltSpliceJunctionFinder.createWriter(mConfig);
                mRetainedIntronWriter = RetainedIntronFinder.createWriter(mConfig);
            }

            if(mConfig.WriteFragmentLengths)
                mFragLengthWriter = FragmentSizeCalcs.createFragmentLengthWriter(mConfig);

            if(mConfig.WriteReadGcRatios)
                mReadGcRatioWriter = GcRatioCounts.createReadGcRatioWriter(mConfig);
        }
    }

    public BufferedWriter getExpRatesWriter() { return mExpRateWriter;}
    public BufferedWriter getAltSpliceJunctionWriter() { return mAltSpliceJunctionWriter;}
    public BufferedWriter getRetainedIntronWriter() { return mRetainedIntronWriter;}
    public BufferedWriter getReadDataWriter() { return mReadDataWriter; }
    public BufferedWriter getFragmentLengthWriter() { return mFragLengthWriter; }
    public BufferedWriter getReadGcRatioWriter() { return mReadGcRatioWriter; }

    public synchronized void writeGeneResult(final GeneResult geneResult)
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
                mGeneDataWriter.write(",TotalFragments,SupportingTrans,Alt,Unspliced,ReadThrough,Chimeric,Duplicates,UnsplicedAlloc,FitResiduals,GeneSet");
                mGeneDataWriter.newLine();
            }

            final EnsemblGeneData geneData = geneResult.geneData();

            long geneLength = geneData.GeneEnd - geneData.GeneStart;

            mGeneDataWriter.write(String.format("%s,%s,%s,%d,%d,%d",
                    geneData.GeneId, geneData.GeneName, geneData.Chromosome, geneLength,
                    geneResult.intronicLength(), geneResult.transCount()));

            mGeneDataWriter.write(String.format(",%d,%d,%d,%d,%d,%d,%d",
                    geneResult.totalFragments(), geneResult.supportingTrans(), geneResult.altFragments(), geneResult.unsplicedFragments(),
                    geneResult.readThroughFragments(), geneResult.chimericFragments(), geneResult.duplicates()));

            mGeneDataWriter.write(String.format(",%.1f,%.1f,%s",
                    geneResult.unsplicedAlloc(), geneResult.fitResiduals(), geneResult.collectionId()));

            mGeneDataWriter.newLine();

        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write gene data file: {}", e.toString());
        }
    }

    public synchronized void writeTranscriptResults(final EnsemblGeneData geneData, final TranscriptResult transResults)
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

            mTransDataWriter.write(String.format("%s,%s,%d,%s,%s,%d,%.0f",
                    geneData.GeneId, geneData.GeneName,
                    transData.TransId, transData.TransName, transData.IsCanonical, transData.exons().size(), transResults.effectiveLength()));

            mTransDataWriter.write(String.format(",%d,%d,%d,%.1f,%d,%d,%.0f",
                    transResults.exonsFound(), transResults.exonicBases(), transResults.exonicBaseCoverage(),
                    transResults.fitAllocation(), transResults.uniqueBases(), transResults.uniqueBaseCoverage(), transResults.uniqueBaseAvgDepth()));

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
            ISF_LOGGER.error("failed to write transcripts data file: {}", e.toString());
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

                final RegionReadData exonReadData = findExonRegion(geneReadData.getExonRegions(), exon.ExonStart, exon.ExonEnd);
                if (exonReadData == null)
                    continue;

                mExonDataWriter.write(String.format("%s,%s,%d,%s",
                        geneReadData.GeneData.GeneId, geneReadData.GeneData.GeneName, transData.TransId, transData.TransName));

                mExonDataWriter.write(String.format(",%d,%d,%d,%d",
                        exon.ExonRank, exon.ExonStart, exon.ExonEnd, exonReadData.getTransExonRefs().size()));

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
            ISF_LOGGER.error("failed to write exon expression file: {}", e.toString());
        }
    }

    public synchronized void writeTransComboCounts(
            final String genesId, final List<String> categories, final double[] counts, final double[] fittedCounts)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mTransComboWriter == null)
            {
                final String outputFileName = mConfig.formOutputFile("category_counts.csv");

                mTransComboWriter = createBufferedWriter(outputFileName, false);
                mTransComboWriter.write("GenesId,Category,Count,FitCount");
                mTransComboWriter.newLine();
            }

            for(int i = 0; i < categories.size(); ++i)
            {
                double count = counts[i];
                final String category = categories.get(i);

                mTransComboWriter.write(String.format("%s,%s,%.0f,%.1f",
                        genesId, category, count, fittedCounts[i]));

                mTransComboWriter.newLine();
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write trans combo data file: {}", e.toString());
        }
    }
}
