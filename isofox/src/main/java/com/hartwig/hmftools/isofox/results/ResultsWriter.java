package com.hartwig.hmftools.isofox.results;

import static com.hartwig.hmftools.common.rna.GeneExpressionFile.GENE_DATA_FILE_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxFunction.NOVEL_LOCATIONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.RETAINED_INTRONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.TRANSCRIPT_COUNTS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.UNMAPPED_READS;
import static com.hartwig.hmftools.isofox.common.FragmentType.ALT;
import static com.hartwig.hmftools.isofox.common.FragmentType.CHIMERIC;
import static com.hartwig.hmftools.isofox.common.FragmentType.DUPLICATE;
import static com.hartwig.hmftools.isofox.common.FragmentType.LOW_MAP_QUAL;
import static com.hartwig.hmftools.isofox.common.FragmentType.TOTAL;
import static com.hartwig.hmftools.isofox.common.FragmentType.TRANS_SUPPORTING;
import static com.hartwig.hmftools.isofox.common.FragmentType.UNSPLICED;
import static com.hartwig.hmftools.isofox.common.FragmentType.typeAsInt;
import static com.hartwig.hmftools.isofox.common.GeneCollection.TRANS_COUNT;
import static com.hartwig.hmftools.isofox.common.GeneCollection.UNIQUE_TRANS_COUNT;
import static com.hartwig.hmftools.isofox.IsofoxFunction.EXPECTED_TRANS_COUNTS;
import static com.hartwig.hmftools.isofox.common.RegionReadData.findExonRegion;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.isofox.BamFragmentAllocator;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.adjusts.FragmentSizeCalcs;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.expression.ExpectedRatesGenerator;
import com.hartwig.hmftools.isofox.expression.TranscriptExpression;
import com.hartwig.hmftools.isofox.adjusts.GcRatioCounts;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunctionFinder;
import com.hartwig.hmftools.isofox.novel.RetainedIntronFinder;
import com.hartwig.hmftools.isofox.novel.SpliceSiteCounter;
import com.hartwig.hmftools.isofox.unmapped.UmrFinder;

public class ResultsWriter
{
    public static final String TRANSCRIPT_RESULTS_FILE = "transcript_data.csv";
    public static final String SPLICE_SITE_FILE = "splice_site_data.csv";
    public static final String SUMMARY_FILE = "summary.csv";

    private final IsofoxConfig mConfig;

    private BufferedWriter mGeneDataWriter;
    private BufferedWriter mGeneCollectionWriter;
    private BufferedWriter mTransDataWriter;
    private BufferedWriter mExonDataWriter;
    private BufferedWriter mCategoryCountsWriter;

    // controlled by other components but instantiated once for output synchronosation
    private BufferedWriter mExpRateWriter;
    private BufferedWriter mReadDataWriter;
    private BufferedWriter mAltSpliceJunctionWriter;
    private BufferedWriter mGeneFragLengthWriter;
    private BufferedWriter mReadGcRatioWriter;
    private BufferedWriter mRetainedIntronWriter;
    private BufferedWriter mSpliceSiteWriter;
    private BufferedWriter mUnamppedReadsWriter;

    public static final String DELIMITER = ",";
    public static final String ITEM_DELIM = ";";
    public static final String SUB_ITEM_DELIM = ":";

    public ResultsWriter(final IsofoxConfig config)
    {
        mConfig = config;

        mGeneDataWriter = null;
        mGeneCollectionWriter = null;
        mTransDataWriter = null;
        mExonDataWriter = null;
        mCategoryCountsWriter = null;
        mExpRateWriter = null;
        mReadDataWriter = null;
        mAltSpliceJunctionWriter = null;
        mGeneFragLengthWriter = null;
        mReadGcRatioWriter = null;
        mRetainedIntronWriter = null;
        mSpliceSiteWriter = null;
        mUnamppedReadsWriter = null;

        if(mConfig.runFunction(TRANSCRIPT_COUNTS))
            initialiseGeneCollectionWriter();

        if(!mConfig.runFusionsOnly())
        {
            initialiseExternalWriters();
        }
    }

    public void close()
    {
        closeBufferedWriter(mGeneDataWriter);
        closeBufferedWriter(mGeneCollectionWriter);
        closeBufferedWriter(mTransDataWriter);
        closeBufferedWriter(mExonDataWriter);
        closeBufferedWriter(mCategoryCountsWriter);
        closeBufferedWriter(mExpRateWriter);
        closeBufferedWriter(mReadDataWriter);
        closeBufferedWriter(mAltSpliceJunctionWriter);
        closeBufferedWriter(mGeneFragLengthWriter);
        closeBufferedWriter(mReadGcRatioWriter);
        closeBufferedWriter(mRetainedIntronWriter);
        closeBufferedWriter(mSpliceSiteWriter);
        closeBufferedWriter(mUnamppedReadsWriter);
    }

    private void initialiseExternalWriters()
    {
        if(mConfig.OutputDir == null)
            return;

        if(mConfig.runFunction(EXPECTED_TRANS_COUNTS) || mConfig.WriteExpectedRates)
        {
            mExpRateWriter = ExpectedRatesGenerator.createWriter(mConfig);
        }

        if(mConfig.WriteFragmentLengthsByGene)
        {
            mGeneFragLengthWriter = FragmentSizeCalcs.createGeneFragmentLengthWriter(mConfig);
            return;
        }

        if(!mConfig.generateExpectedDataOnly())
        {
            if(mConfig.WriteReadData)
                mReadDataWriter = BamFragmentAllocator.createReadDataWriter(mConfig);

            if(mConfig.WriteSpliceSiteData)
                mSpliceSiteWriter = SpliceSiteCounter.createWriter(mConfig);

            if(mConfig.runFunction(NOVEL_LOCATIONS))
                mAltSpliceJunctionWriter = AltSpliceJunctionFinder.createWriter(mConfig);

            if(mConfig.runFunction(RETAINED_INTRONS))
                mRetainedIntronWriter = RetainedIntronFinder.createWriter(mConfig);

            if(mConfig.WriteGcData)
                mReadGcRatioWriter = GcRatioCounts.createReadGcRatioWriter(mConfig);

            if(mConfig.WriteTransComboData)
                mCategoryCountsWriter = TranscriptExpression.createWriter(mConfig);

            if(mConfig.runFunction(UNMAPPED_READS))
                mUnamppedReadsWriter = UmrFinder.createWriter(mConfig);
        }
    }

    public BufferedWriter getExpRatesWriter() { return mExpRateWriter;}
    public BufferedWriter getCategoryCountsWriter() { return mCategoryCountsWriter;}
    public BufferedWriter getAltSpliceJunctionWriter() { return mAltSpliceJunctionWriter;}
    public BufferedWriter getRetainedIntronWriter() { return mRetainedIntronWriter;}
    public BufferedWriter getReadDataWriter() { return mReadDataWriter; }
    public BufferedWriter getSpliceSiteWriter() { return mSpliceSiteWriter; }
    public BufferedWriter getFragmentLengthWriter() { return mGeneFragLengthWriter; }
    public BufferedWriter getReadGcRatioWriter() { return mReadGcRatioWriter; }
    public BufferedWriter getUnmappedReadsWriter() { return mUnamppedReadsWriter; }

    public void writeSummaryStats(final RnaStatistics summaryStats)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            final String outputFileName = mConfig.formOutputFile(SUMMARY_FILE);
            final BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write(RnaStatistics.csvHeader());
            writer.newLine();

            writer.write(summaryStats.toCsv(mConfig.SampleId));
            writer.newLine();
            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write summary data file: {}", e.toString());
        }
    }

    public synchronized void writeGeneResult(final GeneResult geneResult)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mGeneDataWriter == null)
            {
                final String outputFileName = mConfig.formOutputFile(GENE_DATA_FILE_ID);

                mGeneDataWriter = createBufferedWriter(outputFileName, false);
                mGeneDataWriter.write(GeneResult.csvHeader());
                mGeneDataWriter.newLine();
            }

            mGeneDataWriter.write(geneResult.toCsv());
            mGeneDataWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write gene data file: {}", e.toString());
        }
    }

    private void initialiseGeneCollectionWriter()
    {
        if(mConfig.OutputDir == null || mConfig.SampleId == null)
            return;

        try
        {
            final String outputFileName = mConfig.formOutputFile("gene_collection.csv");

            mGeneCollectionWriter = createBufferedWriter(outputFileName, false);
            mGeneCollectionWriter.write("GeneSetId,GeneCount,Chromosome,RangeStart,RangeEnd");
            mGeneCollectionWriter.write(",TotalFragments,Duplicates,SupportingTrans,Unspliced,AltSJ,Chimeric,LowMapQual");
            mGeneCollectionWriter.write(",Genes");
            mGeneCollectionWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write gene collection file: {}", e.toString());
        }
    }

    public synchronized void writeGeneCollectionData(final GeneCollection geneCollection)
    {
        if(mGeneCollectionWriter == null)
            return;

        try
        {
            mGeneCollectionWriter.write(String.format("%s,%d,%s,%d,%d",
                    geneCollection.chrId(), geneCollection.genes().size(), geneCollection.chromosome(),
                    geneCollection.regionBounds()[SE_START], geneCollection.regionBounds()[SE_END]));

            final int[] fragmentCounts = geneCollection.getCounts();

            mGeneCollectionWriter.write(String.format(",%d,%d,%d,%d,%d,%d,%d",
                    fragmentCounts[typeAsInt(TOTAL)], fragmentCounts[typeAsInt(DUPLICATE)], fragmentCounts[typeAsInt(TRANS_SUPPORTING)],
                    fragmentCounts[typeAsInt(UNSPLICED)], fragmentCounts[typeAsInt(ALT)],
                    fragmentCounts[typeAsInt(CHIMERIC)], fragmentCounts[typeAsInt(LOW_MAP_QUAL)]));

            mGeneCollectionWriter.write(String.format(",%s", geneCollection.geneNames(geneCollection.genes().size())));

            mGeneCollectionWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write gene collection file: {}", e.toString());
        }
    }

    public synchronized void writeTranscriptResults(final GeneData geneData, final TranscriptResult transResults)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mTransDataWriter == null)
            {
                final String outputFileName = mConfig.formOutputFile(TRANSCRIPT_RESULTS_FILE);

                mTransDataWriter = createBufferedWriter(outputFileName, false);
                mTransDataWriter.write(TranscriptResult.csvHeader());
                mTransDataWriter.newLine();
            }

            mTransDataWriter.write(transResults.toCsv(geneData));
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

                final RegionReadData exonReadData = findExonRegion(geneReadData.getExonRegions(), exon.Start, exon.End);
                if (exonReadData == null)
                    continue;

                mExonDataWriter.write(String.format("%s,%s,%d,%s",
                        geneReadData.GeneData.GeneId, geneReadData.GeneData.GeneName, transData.TransId, transData.TransName));

                mExonDataWriter.write(String.format(",%d,%d,%d,%d",
                        exon.Rank, exon.Start, exon.End, exonReadData.getTransExonRefs().size()));

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

}
