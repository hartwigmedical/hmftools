package com.hartwig.hmftools.isofox.unmapped;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.seIndex;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.startEndStr;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.ReadRecord.findOverlappingRegions;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_BOUNDARY;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_INTRON;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.common.TransExonRef;

import htsjdk.samtools.SAMRecord;

public class UnmappedReads
{
    private final IsofoxConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;

    private final BufferedWriter mWriter;

    public static final String UNMAPPED_READS_FILE_ID = "unmapped_reads.csv";

    private static final int MIN_SOFT_CLIP_LENGTH = 20;
    private static final int MAX_INTRON_DISTANCE = 5;

    public UnmappedReads(final IsofoxConfig config, final BufferedWriter writer)
    {
        mConfig = config;
        mGeneTransCache = null;

        mWriter = writer;
    }

    public void processReadRecord(final ReadRecord read, final SAMRecord record, final GeneCollection geneCollection)
    {
        // find reads with sufficient soft-clipping and overlapping a splice junction
        if(!read.containsSoftClipping() || read.hasSuppAlignment())
            return;

        int scSide = read.longestSoftClippedEnd();

        int scLength = scSide == SE_START ?
                read.Cigar.getFirstCigarElement().getLength() : read.Cigar.getLastCigarElement().getLength();

        if(scLength < MIN_SOFT_CLIP_LENGTH)
            return;

        final List<RegionReadData> overlappingRegions = findOverlappingRegions(geneCollection.getExonRegions(), read);

        if(overlappingRegions.isEmpty())
            return;

        read.processOverlappingRegions(overlappingRegions);

        StringJoiner transcriptInfo = new StringJoiner(ITEM_DELIM);
        GeneData geneData = null;
        boolean validSpliceFound = false;

        String firstTransName = "";
        int firstExonRank = -1;
        boolean firstIsSpliceAcceptor = false;
        int firstExonBoundary = -1;
        int firstExonBoundaryDistance = -1;

        for(Map.Entry<RegionReadData,RegionMatchType> entry : read.getMappedRegions().entrySet())
        {
            if(entry.getValue() != EXON_INTRON && entry.getValue() != EXON_BOUNDARY)
                continue;

            // capture exon rank and whether a splice acceptor or donor
            RegionReadData region = entry.getKey();
            TransExonRef transExonRef = region.getTransExonRefs().get(0);

            TranscriptData transData = geneCollection.getTranscripts().stream()
                    .filter(x -> x.TransId == transExonRef.TransId).findFirst().orElse(null);

            if(transData == null)
                continue;

            geneData = geneCollection.genes().stream()
                    .filter(x -> x.GeneData.GeneId.equals(transData.GeneId))
                    .map(x -> x.GeneData)
                    .findFirst().orElse(null);

            // work out distance from SC to
            boolean isSpliceAcceptor = (scSide == SE_START) == transData.posStrand();

            int exonBoundaryDistance = 0;
            int exonBoundary = scSide == SE_START ? region.start() : region.end();
            int readBoundary = read.getCoordsBoundary(scSide);

            if(entry.getValue() == EXON_INTRON)
            {
                // soft-clip must be intronic
                if(scSide == SE_START && readBoundary > exonBoundary)
                    continue;
                else if(scSide == SE_END && readBoundary < exonBoundary)
                    continue;

                exonBoundaryDistance = abs(readBoundary - exonBoundary);

                if(exonBoundaryDistance > MAX_INTRON_DISTANCE)
                    continue;
            }

            if(!validSpliceFound)
            {
                validSpliceFound = true;
                firstTransName = transData.TransName;
                firstIsSpliceAcceptor = isSpliceAcceptor;
                firstExonBoundary = exonBoundary;
                firstExonBoundaryDistance = exonBoundaryDistance;
                firstExonRank = transExonRef.ExonRank;
            }

            transcriptInfo.add(String.format("%s:%d:%s:%d:%d",
                    transData.TransName, transExonRef.ExonRank,
                    isSpliceAcceptor ? "acceptor" : "donor", exonBoundary, exonBoundaryDistance));
        }

        if(!validSpliceFound)
            return;

        // determine average base qual within the SC region
        int totalBaseQual = 0;
        int startIndex = scSide == SE_START ? 0 : read.Length - scLength;
        int endIndex = scSide == SE_START ? scLength : read.Length;
        for(int i = startIndex; i < endIndex; ++i)
        {
            totalBaseQual += record.getBaseQualities()[i];
        }

        double avgBaseQual = totalBaseQual / scLength;

        writeReadData(
                mWriter, read, geneData, transcriptInfo.toString(), scSide, scLength, avgBaseQual,
                firstTransName, firstExonRank, firstIsSpliceAcceptor, firstExonBoundary, firstExonBoundaryDistance);
    }

    public static BufferedWriter createWriter(final IsofoxConfig config)
    {
        try
        {
            final String outputFileName = config.formOutputFile(UNMAPPED_READS_FILE_ID);

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("ReadId,Chromosome,PosStart,PosEnd,Cigar,Orientation,SoftClipLength,SoftClipSide,AvgBaseQual");
            writer.write(",GeneId,GeneName,TransName,ExonRank,SpliceType,ExonBoundary,ExonDistance,TranscriptInfo");
            writer.write(",MateMapped,MateChr,MatePosStart,SoftClipBases");
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to create unmapped reads writer: {}", e.toString());
            return null;
        }
    }

    private synchronized static void writeReadData(
            final BufferedWriter writer, final ReadRecord read, final GeneData geneData, final String transcriptInfo,
            int scSide, int scLength, double avgBaseQual, final String transName, int exonRank, boolean isSpliceAcceptor,
            int exonBoundary, int exonBoundaryDistance)
    {
        try
        {
            writer.write(String.format("%s,%s,%d,%d,%s,%d,%d,%s,%.1f",
                    read.Id, read.Chromosome, read.PosStart, read.PosEnd, read.Cigar.toString(), read.orientation(),
                    scLength, startEndStr(scSide), avgBaseQual));

            writer.write(String.format(",%s,%s,%s,%d,%s,%d,%d,%s",
                    geneData.GeneId, geneData.GeneName, transName, exonRank, isSpliceAcceptor ? "acceptor" : "donor",
                    exonBoundary, exonBoundaryDistance, transcriptInfo));

            String scBases = scSide == SE_START ?
                    Nucleotides.reverseStrandBases(read.ReadBases.substring(0, scLength)) : read.ReadBases.substring(read.ReadBases.length() - scLength);

            writer.write(String.format(",%s,%s,%d,%s",
                    !read.isMateUnmapped(), read.mateChromosome(), read.mateStartPosition(), scBases));

            writer.newLine();
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write unmapped read data: {}", e.toString());
        }
    }

}