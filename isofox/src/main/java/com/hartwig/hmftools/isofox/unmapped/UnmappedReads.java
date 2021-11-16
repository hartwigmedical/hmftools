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

    public void processReadRecord(final ReadRecord read, final GeneCollection geneCollection)
    {
        // find reads with sufficient soft-clipping and overlapping a splice junction
        if(!read.containsSoftClipping())
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

            if(entry.getValue() == EXON_INTRON)
            {
                int readBoundary = read.getCoordsBoundary(scSide);
                int exonBoundary = scSide == SE_START ? region.start() : region.end();
                exonBoundaryDistance = abs(readBoundary - exonBoundary);

                if(exonBoundaryDistance > MAX_INTRON_DISTANCE)
                    continue;
            }

            validSpliceFound = true;

            transcriptInfo.add(String.format("%s:%d:%s:%d",
                    transData.TransName, transExonRef.ExonRank,
                    isSpliceAcceptor ? "acceptor" : "donor", exonBoundaryDistance));
        }

        if(validSpliceFound)
        {
            writeReadData(mWriter, read, geneData, transcriptInfo.toString(), scSide, scLength);
        }
    }

    public static BufferedWriter createWriter(final IsofoxConfig config)
    {
        try
        {
            final String outputFileName = config.formOutputFile(UNMAPPED_READS_FILE_ID);

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("ReadId,Chromosome,PosStart,PosEnd,Cigar,SoftClipLength,SoftClipSide");
            writer.write(",GeneId,GeneName,TranscriptInfo");
            writer.write(",MateMapped,MateChr,MatePosStart,SuppData,SoftClipBases");
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
            int scSide, int scLength)
    {
        try
        {
            writer.write(String.format("%s,%s,%d,%d,%s,%d,%s",
                    read.Id, read.Chromosome, read.PosStart, read.PosEnd, read.Cigar.toString(), scLength, startEndStr(scSide)));

            writer.write(String.format(",%s,%s,%s",
                    geneData.GeneId, geneData.GeneName, transcriptInfo));

            String scBases = scSide == SE_START ?
                    Nucleotides.reverseStrandBases(read.ReadBases.substring(0, scLength)) : read.ReadBases.substring(read.ReadBases.length() - scLength);

            writer.write(String.format(",%s,%s,%d,%s,%s",
                    !read.isMateUnmapped(), read.mateChromosome(), read.mateStartPosition(),
                    read.hasSuppAlignment() ? read.getSuppAlignment() : "NONE", scBases));

            writer.newLine();
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write unmapped read data: {}", e.toString());
        }
    }

}