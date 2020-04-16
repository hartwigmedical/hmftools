package com.hartwig.hmftools.isofox.results;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxFunction.NOVEL_LOCATIONS;
import static com.hartwig.hmftools.isofox.common.FragmentType.ALT;
import static com.hartwig.hmftools.isofox.common.FragmentType.CHIMERIC;
import static com.hartwig.hmftools.isofox.common.FragmentType.DUPLICATE;
import static com.hartwig.hmftools.isofox.common.FragmentType.READ_THROUGH;
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

import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.isofox.BamFragmentAllocator;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.adjusts.FragmentSizeCalcs;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator;
import com.hartwig.hmftools.isofox.exp_rates.TranscriptExpression;
import com.hartwig.hmftools.isofox.adjusts.GcRatioCounts;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunctionFinder;
import com.hartwig.hmftools.isofox.novel.RetainedIntronFinder;

public class ResultsWriter
{
    public static final String GENE_RESULTS_FILE = "gene_data.csv";
    public static final String TRANSCRIPT_RESULTS_FILE = "transcript_data.csv";
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

    public static final String DELIMITER = ",";
    public static final String FLD_SAMPLE_ID = "SampleId";
    public static final String FLD_GENE_ID = "GeneId";
    public static final String FLD_GENE_SET_ID = "GeneSetId";
    public static final String FLD_GENE_NAME = "GeneName";
    public static final String FLD_CHROMOSOME = "Chromosome";
    public static final String FLD_TRANS_ID = "TransId";
    public static final String FLD_TRANS_NAME = "TransName";

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

        initialiseGeneCollectionWriter();
        initialiseExternalWriters();
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

            if(mConfig.runFunction(NOVEL_LOCATIONS))
            {
                mAltSpliceJunctionWriter = AltSpliceJunctionFinder.createWriter(mConfig);
                mRetainedIntronWriter = RetainedIntronFinder.createWriter(mConfig);
            }

            if(mConfig.WriteGcData)
                mReadGcRatioWriter = GcRatioCounts.createReadGcRatioWriter(mConfig);

            if(mConfig.WriteTransComboData)
                mCategoryCountsWriter = TranscriptExpression.createWriter(mConfig);
        }
    }

    public BufferedWriter getExpRatesWriter() { return mExpRateWriter;}
    public BufferedWriter getCategoryCountsWriter() { return mCategoryCountsWriter;}
    public BufferedWriter getAltSpliceJunctionWriter() { return mAltSpliceJunctionWriter;}
    public BufferedWriter getRetainedIntronWriter() { return mRetainedIntronWriter;}
    public BufferedWriter getReadDataWriter() { return mReadDataWriter; }
    public BufferedWriter getFragmentLengthWriter() { return mGeneFragLengthWriter; }
    public BufferedWriter getReadGcRatioWriter() { return mReadGcRatioWriter; }

    public void writeSummaryStats(final SummaryStats summaryStats)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            final String outputFileName = mConfig.formOutputFile(SUMMARY_FILE);
            final BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write(SummaryStats.csvHeader());
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
                final String outputFileName = mConfig.formOutputFile(GENE_RESULTS_FILE);

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
        if(mConfig.OutputDir == null)
            return;

        try
        {
            final String outputFileName = mConfig.formOutputFile("gene_collection_data.csv");

            mGeneCollectionWriter = createBufferedWriter(outputFileName, false);
            mGeneCollectionWriter.write("GeneSetId,GeneCount,Chromosome,RangeStart,RangeEnd");
            mGeneCollectionWriter.write(",TotalFragments,Duplicates,SupportingTrans,Unspliced,AltSJ,ReadThrough,Chimeric");
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
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            mGeneCollectionWriter.write(String.format("%s,%d,%s,%d,%d",
                    geneCollection.chrId(), geneCollection.genes().size(), geneCollection.chromosome(),
                    geneCollection.regionBounds()[SE_START], geneCollection.regionBounds()[SE_END]));

            final int[] fragmentCounts = geneCollection.getCounts();

            mGeneCollectionWriter.write(String.format(",%d,%d,%d,%d,%d,%d,%d",
                    fragmentCounts[typeAsInt(TOTAL)], fragmentCounts[typeAsInt(DUPLICATE)], fragmentCounts[typeAsInt(TRANS_SUPPORTING)],
                    fragmentCounts[typeAsInt(UNSPLICED)], fragmentCounts[typeAsInt(ALT)], fragmentCounts[typeAsInt(READ_THROUGH)],
                    fragmentCounts[typeAsInt(CHIMERIC)]));

            mGeneCollectionWriter.write(String.format(",%s", geneCollection.geneNames(geneCollection.genes().size())));

            mGeneCollectionWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write gene collection file: {}", e.toString());
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

}
