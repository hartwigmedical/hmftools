package com.hartwig.hmftools.isofox.results;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.rna.CanonicalSpliceJunctionFile.CANONICAL_SJ_FILE_ID;
import static com.hartwig.hmftools.common.rna.GeneExpressionFile.GENE_EXPRESSION_FILE_ID;
import static com.hartwig.hmftools.common.rna.GeneExpressionFile.TRANSCRIPT_EXPRESSION_FILE_ID;
import static com.hartwig.hmftools.common.rna.RnaStatistics.SUMMARY_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxFunction.ALT_SPLICE_JUNCTIONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.READ_COUNTS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.RETAINED_INTRONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.TRANSCRIPT_COUNTS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.UNMAPPED_READS;
import static com.hartwig.hmftools.isofox.common.FragmentType.ALT;
import static com.hartwig.hmftools.isofox.common.FragmentType.CHIMERIC;
import static com.hartwig.hmftools.isofox.common.FragmentType.DUPLICATE;
import static com.hartwig.hmftools.isofox.common.FragmentType.FORWARD_STRAND;
import static com.hartwig.hmftools.isofox.common.FragmentType.LOW_MAP_QUAL;
import static com.hartwig.hmftools.isofox.common.FragmentType.REVERSE_STRAND;
import static com.hartwig.hmftools.isofox.common.FragmentType.TOTAL;
import static com.hartwig.hmftools.isofox.common.FragmentType.TRANS_SUPPORTING;
import static com.hartwig.hmftools.isofox.common.FragmentType.UNSPLICED;
import static com.hartwig.hmftools.isofox.common.GeneCollection.TRANS_COUNT;
import static com.hartwig.hmftools.isofox.common.GeneCollection.UNIQUE_TRANS_COUNT;
import static com.hartwig.hmftools.isofox.common.RegionReadData.findExonRegion;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.rna.CanonicalSpliceJunctionFile;
import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.isofox.FragmentAllocator;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.adjusts.FragmentSizeCalcs;
import com.hartwig.hmftools.isofox.common.BamReadCounter;
import com.hartwig.hmftools.isofox.common.FragmentTypeCounts;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.common.TransExonRef;
import com.hartwig.hmftools.isofox.expression.ExpectedRatesCommon;
import com.hartwig.hmftools.isofox.expression.TranscriptExpression;
import com.hartwig.hmftools.isofox.adjusts.GcRatioCounts;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunctionFinder;
import com.hartwig.hmftools.isofox.novel.RetainedIntronFinder;
import com.hartwig.hmftools.isofox.novel.SpliceSiteCounter;
import com.hartwig.hmftools.isofox.unmapped.UmrFinder;

import com.google.common.collect.Lists;

public class ResultsWriter
{
    public static final String SPLICE_SITE_FILE = "splice_site_data.csv";

    private final IsofoxConfig mConfig;

    private BufferedWriter mGeneDataWriter;
    private BufferedWriter mGeneCollectionWriter;
    private BufferedWriter mTransDataWriter;
    private BufferedWriter mExonDataWriter;
    private BufferedWriter mSpliceJunctionWriter;
    private BufferedWriter mCategoryCountsWriter;

    // controlled by other components but instantiated once for output synchronosation
    private BufferedWriter mReadDataWriter;
    private BufferedWriter mAltSpliceJunctionWriter;
    private BufferedWriter mGeneFragLengthWriter;
    private BufferedWriter mReadGcRatioWriter;
    private BufferedWriter mRetainedIntronWriter;
    private BufferedWriter mSpliceSiteWriter;
    private BufferedWriter mUnamppedReadsWriter;

    public static final String DELIMITER = CSV_DELIM;

    public ResultsWriter(final IsofoxConfig config)
    {
        mConfig = config;

        mGeneDataWriter = null;
        mGeneCollectionWriter = null;
        mTransDataWriter = null;
        mExonDataWriter = null;
        mSpliceJunctionWriter = null;
        mCategoryCountsWriter = null;
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
        closeBufferedWriter(mSpliceJunctionWriter);
        closeBufferedWriter(mCategoryCountsWriter);
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

        if(mConfig.WriteFragmentLengthsByGene)
        {
            mGeneFragLengthWriter = FragmentSizeCalcs.createGeneFragmentLengthWriter(mConfig);
            return;
        }

        if(mConfig.WriteReadData)
        {
            if(mConfig.runFunction(READ_COUNTS))
                mReadDataWriter = BamReadCounter.createReadDataWriter(mConfig);
            else
                mReadDataWriter = FragmentAllocator.createReadDataWriter(mConfig);
        }

        if(mConfig.WriteSpliceSiteData)
            mSpliceSiteWriter = SpliceSiteCounter.createWriter(mConfig);

        if(mConfig.runFunction(ALT_SPLICE_JUNCTIONS))
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
            final String outputFileName = mConfig.formOutputFile(SUMMARY_FILE_ID);
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
                final String outputFileName = mConfig.formOutputFile(GENE_EXPRESSION_FILE_ID);

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
            mGeneCollectionWriter.write(",ForwardStrand,ReverseStrand");
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
            mGeneCollectionWriter.write(format("%s,%d,%s,%d,%d",
                    geneCollection.chrId(), geneCollection.genes().size(), geneCollection.chromosome(),
                    geneCollection.regionBounds()[SE_START], geneCollection.regionBounds()[SE_END]));

            final FragmentTypeCounts fragmentCounts = geneCollection.fragmentTypeCounts();

            mGeneCollectionWriter.write(format(",%d,%d,%d,%d,%d,%d,%d,%d,%d",
                    fragmentCounts.typeCount(TOTAL), fragmentCounts.typeCount(DUPLICATE), fragmentCounts.typeCount(TRANS_SUPPORTING),
                    fragmentCounts.typeCount(UNSPLICED), fragmentCounts.typeCount(ALT),
                    fragmentCounts.typeCount(CHIMERIC), fragmentCounts.typeCount(LOW_MAP_QUAL),
                    fragmentCounts.typeCount(FORWARD_STRAND), fragmentCounts.typeCount(REVERSE_STRAND)));

            mGeneCollectionWriter.write(format(",%s", geneCollection.geneNames(geneCollection.genes().size())));

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
                final String outputFileName = mConfig.formOutputFile(TRANSCRIPT_EXPRESSION_FILE_ID);

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
                if(exonReadData == null)
                    continue;

                mExonDataWriter.write(format("%s,%s,%d,%s",
                        geneReadData.GeneData.GeneId, geneReadData.GeneData.GeneName, transData.TransId, transData.TransName));

                mExonDataWriter.write(format(",%d,%d,%d,%d",
                        exon.Rank, exon.Start, exon.End, exonReadData.getTransExonRefs().size()));

                int[] matchCounts = exonReadData.getTranscriptReadCount(transData.TransId);
                int[] startSjCounts = exonReadData.getTranscriptJunctionMatchCount(transData.TransId, SE_START);
                int[] endSjCounts = exonReadData.getTranscriptJunctionMatchCount(transData.TransId, SE_END);

                int uniqueBaseTotalDepth = exonReadData.uniqueBaseTotalDepth();
                int uniqueBaseCount = exonReadData.uniqueBaseCount();
                double uniqueAvgDepth = uniqueBaseCount > 0 ? uniqueBaseTotalDepth / (double)uniqueBaseCount : 0;

                mExonDataWriter.write(format(",%d,%.0f,%d,%d,%.0f",
                        exonReadData.baseCoverage(1), exonReadData.averageDepth(),
                        uniqueBaseCount, exonReadData.uniqueBaseCoverage(1), uniqueAvgDepth));

                mExonDataWriter.write(format(",%d,%d,%d,%d,%d,%d",
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

    private static final int MIN_SPLICE_JUNCTON_FRAGMENTS = 3;

    public synchronized void writeSpliceJunctionData(final GeneCollection geneCollection)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mSpliceJunctionWriter == null)
            {
                final String outputFileName = mConfig.formOutputFile(CANONICAL_SJ_FILE_ID);

                mSpliceJunctionWriter = createBufferedWriter(outputFileName, false);
                mSpliceJunctionWriter.write(CanonicalSpliceJunctionFile.csvHeader());
                mSpliceJunctionWriter.newLine();
            }

            Map<String,SpliceJunctionData> junctionCounts = Maps.newHashMap();

            for(RegionReadData regionReadData : geneCollection.getExonRegions())
            {
                if(regionReadData.getPreRegions().isEmpty())
                    continue;

                int nextExonStart = regionReadData.start();

                Map<Integer,int[][]> nextTransJunctionCounts = regionReadData.getTranscriptJunctionCounts();

                for(RegionReadData prevRegion : regionReadData.getPreRegions())
                {
                    int prevExonEnd = prevRegion.end();

                    String junctionStr = format("%d_%d", prevExonEnd, nextExonStart);

                    if(junctionCounts.containsKey(junctionStr))
                        continue;

                    Map<Integer, int[][]> prevTransJunctionCounts = prevRegion.getTranscriptJunctionCounts();

                    SpliceJunctionData sjData = null;

                    for(Map.Entry<Integer, int[][]> entry : nextTransJunctionCounts.entrySet())
                    {
                        Integer nextTransId = entry.getKey();

                        final int[][] prevCounts = prevTransJunctionCounts.get(nextTransId);

                        if(prevCounts == null)
                            continue;

                        TransExonRef transExonRef = regionReadData.getTransExonRefs().stream()
                                .filter(x -> x.TransId == nextTransId).findFirst().orElse(null);

                        if(transExonRef == null)
                            continue;

                        if(sjData == null)
                        {
                            final int[][] nextCounts = entry.getValue();
                            int prevSpliceCount = prevCounts[SE_END][TRANS_COUNT];
                            int nextSpliceCount = nextCounts[SE_START][TRANS_COUNT];

                            if(prevSpliceCount != nextSpliceCount) // they should have been incremented equally
                                continue;

                            GeneReadData geneData = geneCollection.genes().stream()
                                    .filter(x -> x.GeneData.GeneId.equals(transExonRef.GeneId)).findFirst().orElse(null);

                            sjData = new SpliceJunctionData(
                                    geneData.GeneData.GeneId, geneData.GeneData.GeneName,
                                    prevExonEnd, nextExonStart, prevRegion.getBoundaryBaseDepth(SE_END),
                                    regionReadData.getBoundaryBaseDepth(SE_START), nextSpliceCount);

                            junctionCounts.put(junctionStr, sjData);
                        }

                        sjData.TranscriptNames.add(transExonRef.TransName);
                    }
                }
            }

            for(Map.Entry<String,SpliceJunctionData> entry : junctionCounts.entrySet())
            {
                SpliceJunctionData sjData = entry.getValue();

                if(sjData.FragmentCount < MIN_SPLICE_JUNCTON_FRAGMENTS)
                    continue;

                StringJoiner sj = new StringJoiner(ITEM_DELIM);
                sjData.TranscriptNames.forEach(x -> sj.add(x));

                mSpliceJunctionWriter.write(format("%s,%s,%s,%s,%s,%d,%d,%d,%s",
                        sjData.GeneId, sjData.GeneName, geneCollection.chromosome(),
                        sjData.SjStart, sjData.SjEnd, sjData.FragmentCount, sjData.DepthStart, sjData.DepthEnd, sj.toString()));

                mSpliceJunctionWriter.newLine();
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write splice junction file: {}", e.toString());
        }
    }

    private class SpliceJunctionData
    {
        public final String GeneId;
        public final String GeneName;
        public final int SjStart;
        public final int SjEnd;
        public final int DepthStart;
        public final int DepthEnd;
        public final int FragmentCount;
        public final List<String> TranscriptNames;

        public SpliceJunctionData(
                final String geneId, final String geneName,
                final int sjStart, final int sjEnd, final int depthStart, final int depthEnd, final int fragmentCount)
        {
            GeneId = geneId;
            GeneName = geneName;
            SjStart = sjStart;
            SjEnd = sjEnd;
            DepthStart = depthStart;
            DepthEnd = depthEnd;
            FragmentCount = fragmentCount;
            TranscriptNames = Lists.newArrayList();
        }
    }
}
