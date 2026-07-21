package com.hartwig.hmftools.isofox.results;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SupplementaryReadData.fromAlignment;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.isofox.WriteType.CHIMERIC_POSITION_DATA;
import static com.hartwig.hmftools.isofox.WriteType.CHIMERIC_READ;
import static com.hartwig.hmftools.isofox.WriteType.FRAG_LENGTH_BY_GENE;
import static com.hartwig.hmftools.isofox.WriteType.GC_RATIO;
import static com.hartwig.hmftools.isofox.WriteType.MULTI_MAP_LOCI;
import static com.hartwig.hmftools.isofox.WriteType.READ;
import static com.hartwig.hmftools.isofox.WriteType.SPLICE_SITE;
import static com.hartwig.hmftools.isofox.WriteType.TRANS_COMBO;
import static com.hartwig.hmftools.isofox.common.Read.clippedSide;
import static com.hartwig.hmftools.isofox.novel.CanonicalSpliceJunctionFile.CANONICAL_SJ_FILE_ID;
import static com.hartwig.hmftools.common.rna.GeneExpressionFile.GENE_EXPRESSION_FILE_ID;
import static com.hartwig.hmftools.common.rna.RnaStatisticFile.SUMMARY_FILE_ID;
import static com.hartwig.hmftools.common.rna.TranscriptExpressionFile.TRANSCRIPT_EXPRESSION_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxFunction.ALT_SPLICE_JUNCTIONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.READ_COUNTS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.RETAINED_INTRONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.TRANSCRIPT_COUNTS;
import static com.hartwig.hmftools.isofox.common.FragmentType.ALT;
import static com.hartwig.hmftools.isofox.common.FragmentType.CHIMERIC;
import static com.hartwig.hmftools.isofox.common.FragmentType.DUPLICATE;
import static com.hartwig.hmftools.isofox.common.FragmentType.FORWARD_STRAND;
import static com.hartwig.hmftools.isofox.common.FragmentType.MULTI_MAPPED;
import static com.hartwig.hmftools.isofox.common.FragmentType.REVERSE_STRAND;
import static com.hartwig.hmftools.isofox.common.FragmentType.TOTAL;
import static com.hartwig.hmftools.isofox.common.FragmentType.TRANS_SUPPORTING;
import static com.hartwig.hmftools.isofox.common.FragmentType.UNSPLICED;
import static com.hartwig.hmftools.isofox.common.GeneCollection.TRANS_COUNT;
import static com.hartwig.hmftools.isofox.common.GeneCollection.UNIQUE_TRANS_COUNT;
import static com.hartwig.hmftools.isofox.common.RegionReadData.findExonRegion;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.ClippedSide;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.isofox.WriteType;
import com.hartwig.hmftools.isofox.common.BaseDepth;
import com.hartwig.hmftools.isofox.common.Read;
import com.hartwig.hmftools.isofox.fusion.ChimericPosData;
import com.hartwig.hmftools.isofox.fusion.ChimericReadGroup;
import com.hartwig.hmftools.isofox.fusion.ChimericRemoteRegion;
import com.hartwig.hmftools.isofox.novel.CanonicalSpliceJunctionFile;
import com.hartwig.hmftools.common.rna.RnaStatisticFile;
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
import com.hartwig.hmftools.isofox.expression.TranscriptExpression;
import com.hartwig.hmftools.isofox.adjusts.GcRatioCounts;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunctionFinder;
import com.hartwig.hmftools.isofox.novel.RetainedIntronFinder;
import com.hartwig.hmftools.isofox.novel.SpliceSiteCounter;

import com.google.common.collect.Lists;

public class ResultsWriter
{
    public static final String SPLICE_SITE_FILE = "splice_site_data.csv";

    public static final String OLD_FILE_DELIM = CSV_DELIM; // for cohort files and non-standard output for now

    private final IsofoxConfig mConfig;

    private BufferedWriter mGeneDataWriter;
    private BufferedWriter mGeneCollectionWriter;
    private BufferedWriter mTransDataWriter;
    private BufferedWriter mExonDataWriter;
    private BufferedWriter mSpliceJunctionWriter;
    private BufferedWriter mCategoryCountsWriter;

    // controlled by other components but instantiated once for output synchronisation
    private BufferedWriter mReadDataWriter;
    private BufferedWriter mAltSjUnfilteredWriter;
    private BufferedWriter mAltSjPassingWriter;
    private BufferedWriter mGeneFragLengthWriter;
    private BufferedWriter mReadGcRatioWriter;
    private BufferedWriter mRetainedIntronWriter;
    private BufferedWriter mSpliceSiteWriter;
    private BufferedWriter mChimericReadWriter;
    private BufferedWriter mChimericPositionDataWriter;
    private BufferedWriter mMultiMapLociWriter;

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
        mAltSjUnfilteredWriter = null;
        mAltSjPassingWriter = null;
        mGeneFragLengthWriter = null;
        mReadGcRatioWriter = null;
        mRetainedIntronWriter = null;
        mSpliceSiteWriter = null;
        mChimericReadWriter = null;
        mChimericPositionDataWriter = null;

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
        closeBufferedWriter(mAltSjUnfilteredWriter);
        closeBufferedWriter(mAltSjPassingWriter);
        closeBufferedWriter(mGeneFragLengthWriter);
        closeBufferedWriter(mReadGcRatioWriter);
        closeBufferedWriter(mRetainedIntronWriter);
        closeBufferedWriter(mSpliceSiteWriter);
        closeBufferedWriter(mChimericReadWriter);
        closeBufferedWriter(mChimericPositionDataWriter);
        closeBufferedWriter(mMultiMapLociWriter);
    }

    private void initialiseExternalWriters()
    {
        if(mConfig.OutputDir == null)
            return;

        if(mConfig.writeType(FRAG_LENGTH_BY_GENE))
        {
            mGeneFragLengthWriter = FragmentSizeCalcs.createGeneFragmentLengthWriter(mConfig);
            return;
        }

        if(mConfig.writeType(READ))
        {
            if(mConfig.runFunction(READ_COUNTS))
                mReadDataWriter = BamReadCounter.createReadDataWriter(mConfig);
            else
                mReadDataWriter = FragmentAllocator.createReadDataWriter(mConfig);
        }

        if(mConfig.writeType(SPLICE_SITE))
            mSpliceSiteWriter = SpliceSiteCounter.createWriter(mConfig);

        if(mConfig.runFunction(ALT_SPLICE_JUNCTIONS))
        {
            mAltSjUnfilteredWriter = AltSpliceJunctionFinder.createUnfilteredWriter(mConfig);
            mAltSjPassingWriter = AltSpliceJunctionFinder.createPassingWriter(mConfig);
        }

        if(mConfig.runFunction(RETAINED_INTRONS))
            mRetainedIntronWriter = RetainedIntronFinder.createWriter(mConfig);

        if(mConfig.writeType(GC_RATIO))
            mReadGcRatioWriter = GcRatioCounts.createReadGcRatioWriter(mConfig);

        if(mConfig.writeType(TRANS_COMBO))
            mCategoryCountsWriter = TranscriptExpression.createWriter(mConfig);

        if(mConfig.writeType(CHIMERIC_READ))
            initialiseChimericReadWriter();

        if(mConfig.writeType(CHIMERIC_POSITION_DATA))
            initialiseChimericPositionDataWriter();

        if(mConfig.writeType(MULTI_MAP_LOCI))
            mMultiMapLociWriter = FragmentAllocator.createMultiMapLociWriter(mConfig);
    }

    public BufferedWriter getCategoryCountsWriter() { return mCategoryCountsWriter;}
    public BufferedWriter getAltSjUnfilteredWriter() { return mAltSjUnfilteredWriter;}
    public BufferedWriter getAltSjPassingWriter() { return mAltSjPassingWriter;}
    public BufferedWriter getRetainedIntronWriter() { return mRetainedIntronWriter;}
    public BufferedWriter getReadDataWriter() { return mReadDataWriter; }
    public BufferedWriter getMultiMapLociWriter() { return mMultiMapLociWriter; }
    public BufferedWriter getSpliceSiteWriter() { return mSpliceSiteWriter; }
    public BufferedWriter getFragmentLengthWriter() { return mGeneFragLengthWriter; }
    public BufferedWriter getReadGcRatioWriter() { return mReadGcRatioWriter; }
    public BufferedWriter getChimericReadWriter() { return mChimericReadWriter; }
    public BufferedWriter getChimericPositionDataWriter() { return mChimericPositionDataWriter; }

    public void writeSummaryStats(final RnaStatistics summaryStats)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            String outputFileName = mConfig.formOutputFile(SUMMARY_FILE_ID);
            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write(RnaStatisticFile.header());
            writer.newLine();

            writer.write(RnaStatisticFile.writeLine(mConfig.SampleId, summaryStats));
            writer.newLine();
            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write summary data file: {}", e.toString());
        }
    }

    public synchronized void writeGeneExpression(final GeneResult geneResult)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mGeneDataWriter == null)
            {
                String outputFileName = mConfig.formOutputFile(GENE_EXPRESSION_FILE_ID);

                mGeneDataWriter = createBufferedWriter(outputFileName, false);
                mGeneDataWriter.write(GeneResult.header());
                mGeneDataWriter.newLine();
            }

            mGeneDataWriter.write(geneResult.toLine());
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
            final String outputFileName = mConfig.formOutputFile("gene_collection.tsv");

            mGeneCollectionWriter = createBufferedWriter(outputFileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("GeneSetId").add("GeneCount").add("Chromosome").add("RangeStart").add("RangeEnd");
            sj.add("TotalFragments").add("Duplicates").add("SupportingTrans").add("Unspliced").add("AltSJ").add("Chimeric");
            sj.add("LowMapQual").add("ForwardStrand").add("ReverseStrand").add("Genes");

            mGeneCollectionWriter.write(sj.toString());
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
            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(geneCollection.chrId());
            sj.add(String.valueOf(geneCollection.genes().size()));
            sj.add(geneCollection.chromosome());
            sj.add(String.valueOf(geneCollection.regionBounds()[SE_START]));
            sj.add(String.valueOf(geneCollection.regionBounds()[SE_END]));

            final FragmentTypeCounts fragmentCounts = geneCollection.fragmentTypeCounts();
            sj.add(String.valueOf(fragmentCounts.typeCount(TOTAL)));
            sj.add(String.valueOf(fragmentCounts.typeCount(DUPLICATE)));
            sj.add(String.valueOf(fragmentCounts.typeCount(TRANS_SUPPORTING)));
            sj.add(String.valueOf(fragmentCounts.typeCount(UNSPLICED)));
            sj.add(String.valueOf(fragmentCounts.typeCount(ALT)));
            sj.add(String.valueOf(fragmentCounts.typeCount(CHIMERIC)));
            sj.add(String.valueOf(fragmentCounts.typeCount(MULTI_MAPPED)));
            sj.add(String.valueOf(fragmentCounts.typeCount(FORWARD_STRAND)));
            sj.add(String.valueOf(fragmentCounts.typeCount(REVERSE_STRAND)));

            sj.add(geneCollection.geneNames(geneCollection.genes().size()));

            mGeneCollectionWriter.write(sj.toString());

            mGeneCollectionWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write gene collection file: {}", e.toString());
        }
    }

    public synchronized void writeTranscriptExpression(final GeneData geneData, final TranscriptResult transResults)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mTransDataWriter == null)
            {
                final String outputFileName = mConfig.formOutputFile(TRANSCRIPT_EXPRESSION_FILE_ID);

                mTransDataWriter = createBufferedWriter(outputFileName, false);
                mTransDataWriter.write(TranscriptResult.header());
                mTransDataWriter.newLine();
            }

            mTransDataWriter.write(transResults.toLine(geneData));
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
                        geneReadData.Gene.GeneId, geneReadData.Gene.GeneName, transData.TransId, transData.TransName));

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

    private void initialiseChimericReadWriter()
    {
        try
        {
            final String outputFileName = mConfig.formOutputFile("chimeric_reads.tsv");
            mChimericReadWriter = createBufferedWriter(outputFileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add("GroupCount").add("GroupComplete").add("ReadId");
            sj.add("Chromosome").add("PosStart").add("PosEnd").add("Cigar").add("Flags").add("MapQual");
            sj.add("IsSupp").add("IsDup").add("MateChr").add("MatePosition").add("SuppChr").add("SuppPosition");
            sj.add("GeneSet").add("GeneName").add("BaseDepth");

            mChimericReadWriter.write(sj.toString());
            mChimericReadWriter.newLine();
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to initialise chimeric read data: {}", e.toString());
        }
    }

    public static synchronized void writeChimericReadData(
            final BufferedWriter writer, final ChimericReadGroup readGroup, final BaseDepth baseDepth)
    {
        if(writer == null)
            return;

        try
        {
            Read primaryRead = null;
            ClippedSide maxClippedSide = null;
            String suppChromosome = "";
            int suppPosition = 0;

            for(Read read : readGroup.reads())
            {
                ClippedSide clippedSide = clippedSide(read);

                if(primaryRead == null)
                {
                    primaryRead = read;
                    maxClippedSide = clippedSide;
                }
                else if(primaryRead.isSupplementaryAlignment() && !read.isSupplementaryAlignment())
                {
                    primaryRead = read;
                    maxClippedSide = clippedSide;
                }
                else
                {

                    if(clippedSide.Length > maxClippedSide.Length)
                    {
                        primaryRead = read;
                        maxClippedSide = clippedSide;
                    }
                }

                if(suppChromosome.isEmpty() && !read.isSupplementaryAlignment() && read.hasSuppAlignment())
                {
                    String[] suppDataItems = read.getSuppAlignment().split(CSV_DELIM, -1);

                    if(suppDataItems.length > 2)
                    {
                        suppChromosome = suppDataItems[0];
                        suppPosition = Integer.parseInt(suppDataItems[1]);
                    }
                }
            }

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(String.valueOf(readGroup.size()));
            sj.add(String.valueOf(readGroup.isComplete()));
            sj.add(primaryRead.Id);
            sj.add(primaryRead.Chromosome);
            sj.add(String.valueOf(primaryRead.PosStart));
            sj.add(String.valueOf(primaryRead.PosEnd));
            sj.add(primaryRead.cigarStr());
            sj.add(String.valueOf(primaryRead.flags()));
            sj.add(String.valueOf(primaryRead.mapQuality()));
            sj.add(String.valueOf(primaryRead.isSupplementaryAlignment()));
            sj.add(String.valueOf(primaryRead.isDuplicate()));
            sj.add(primaryRead.mateChromosome());
            sj.add(String.valueOf(primaryRead.mateStartPosition()));

            sj.add(suppChromosome);
            sj.add(String.valueOf(suppPosition));

            sj.add(String.valueOf(primaryRead.getGeneCollectons()[SE_START]));

            String geneId = "";

            if(!primaryRead.getReadTransExonRefs().isEmpty())
            {
                TransExonRef transExonRef = primaryRead.getReadTransExonRefs().values().iterator().next().get(0);
                geneId = transExonRef.GeneId;
            }

            sj.add(geneId);

            int basePosition = maxClippedSide.Length > 0 ?
                    (maxClippedSide.isLeft() ? primaryRead.PosStart : primaryRead.PosEnd) : primaryRead.PosStart;

            sj.add(String.valueOf(baseDepth.depthAtBase(basePosition)));

            writer.write(sj.toString());
            writer.newLine();
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write chimeric read data: {}", e.toString());
        }
    }

    private void initialiseChimericPositionDataWriter()
    {
        try
        {
            final String outputFileName = mConfig.formOutputFile("chimeric_pos_data.tsv");
            mChimericPositionDataWriter = createBufferedWriter(outputFileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add("Chromosome").add("Position").add("FragCount").add("DupCount").add("SuppCount");
            sj.add("RemoteMaxRegion").add("RemoteMaxCount").add("RemoteTotalRegions");
            sj.add("InitialReadId").add("InitialReadBases");

            mChimericPositionDataWriter.write(sj.toString());
            mChimericPositionDataWriter.newLine();
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to initialise chimeric position data: {}", e.toString());
        }
    }

    private static final int MIN_CHRIMERIC_POS_FRAG_COUNT = 5;

    public static synchronized void writeChimericPositionData(final BufferedWriter writer, final Collection<ChimericPosData> posDataList)
    {
        if(writer == null)
            return;

        try
        {
            for(ChimericPosData posData : posDataList)
            {
                if(posData.ReadCount < MIN_CHRIMERIC_POS_FRAG_COUNT)
                    continue;

                StringJoiner sj = new StringJoiner(TSV_DELIM);
                sj.add(posData.Chromosome);
                sj.add(String.valueOf(posData.Position));
                sj.add(String.valueOf(posData.ReadCount));
                sj.add(String.valueOf(posData.DuplicateCount));
                sj.add(String.valueOf(posData.SuppCount));

                if(!posData.RemoteRegions.isEmpty())
                {
                    Collections.sort(posData.RemoteRegions, Comparator.comparingInt(x -> -x.Count));
                    ChimericRemoteRegion maxRemoteRegion = posData.RemoteRegions.get(0);

                    sj.add(maxRemoteRegion.toString());
                    sj.add(String.valueOf(maxRemoteRegion.Count));
                    sj.add(String.valueOf(posData.RemoteRegions.size()));
                }
                else
                {
                    sj.add("").add("0").add("0");
                }

                sj.add(posData.InitialReadId);
                sj.add(posData.InitialReadBases);

                writer.write(sj.toString());
                writer.newLine();
            }
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write chimeric position data: {}", e.toString());
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
                mSpliceJunctionWriter.write(CanonicalSpliceJunctionFile.header());
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
                                    .filter(x -> x.Gene.GeneId.equals(transExonRef.GeneId)).findFirst().orElse(null);

                            sjData = new SpliceJunctionData(
                                    geneData.Gene.GeneId, geneData.Gene.GeneName,
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
