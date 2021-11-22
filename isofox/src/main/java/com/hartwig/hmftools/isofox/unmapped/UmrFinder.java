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
import static com.hartwig.hmftools.isofox.unmapped.UnmappedRead.SPLICE_TYPE_ACCEPTOR;
import static com.hartwig.hmftools.isofox.unmapped.UnmappedRead.SPLICE_TYPE_DONOR;
import static com.hartwig.hmftools.isofox.unmapped.UnmappedRead.positionKey;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.common.TransExonRef;

import org.apache.commons.compress.utils.Lists;

import htsjdk.samtools.SAMRecord;

public class UmrFinder
{
    private final IsofoxConfig mConfig;
    private final UmrCohortFrequency mCohortFrequency;
    private GeneCollection mGenes;
    private final List<UnmappedRead> mCandidateReads;
    private final Set<String> mSupplementaryReadKeys;

    private final BufferedWriter mWriter;

    public static final String UNMAPPED_READS_FILE_ID = "unmapped_reads.csv";

    private static final int MIN_SOFT_CLIP_LENGTH = 20;
    private static final int MAX_INTRON_DISTANCE = 5;

    public UmrFinder(final IsofoxConfig config, final BufferedWriter writer)
    {
        mConfig = config;
        mCohortFrequency = new UmrCohortFrequency(config.UnmappedCohortFreqFile);

        mCandidateReads = Lists.newArrayList();
        mSupplementaryReadKeys = Sets.newHashSet();
        mGenes = null;

        mWriter = writer;
    }

    public void setGeneData(final GeneCollection genes)
    {
        mGenes = genes;
        mCandidateReads.clear();
        mSupplementaryReadKeys.clear();
    }

    public void processReadRecord(final ReadRecord read, final SAMRecord record)
    {
        // find reads with sufficient soft-clipping and overlapping a splice junction
        if(!read.containsSoftClipping())
            return;

        int scSide = read.longestSoftClippedEnd();

        int scLength = scSide == SE_START ?
                read.Cigar.getFirstCigarElement().getLength() : read.Cigar.getLastCigarElement().getLength();

        if(scLength < MIN_SOFT_CLIP_LENGTH)
            return;

        final List<RegionReadData> overlappingRegions = findOverlappingRegions(mGenes.getExonRegions(), read);

        if(overlappingRegions.isEmpty())
            return;

        read.processOverlappingRegions(overlappingRegions);

        // StringJoiner transcriptInfo = new StringJoiner(ITEM_DELIM);

        TranscriptData selectedTransData = null;
        int exonRank = -1;
        String spliceType = "";
        int exonBoundary = -1;
        int exonBoundaryDistance = -1;

        for(Map.Entry<RegionReadData,RegionMatchType> entry : read.getMappedRegions().entrySet())
        {
            if(entry.getValue() != EXON_INTRON && entry.getValue() != EXON_BOUNDARY)
                continue;

            // capture exon rank and whether a splice acceptor or donor
            RegionReadData region = entry.getKey();
            TransExonRef transExonRef = region.getTransExonRefs().get(0);

            TranscriptData transData = mGenes.getTranscripts().stream()
                    .filter(x -> x.TransId == transExonRef.TransId).findFirst().orElse(null);

            if(transData == null)
                continue;

            // work out distance from SC to
            boolean isSpliceAcceptor = (scSide == SE_START) == transData.posStrand();

            // ignore
            if(isSpliceAcceptor && transExonRef.ExonRank == 1)
                continue;

            if(!isSpliceAcceptor && transExonRef.ExonRank == transData.exons().size())
                continue;

            int boundaryDistance = 0;
            int regionBoundary = scSide == SE_START ? region.start() : region.end();
            int readBoundary = read.getCoordsBoundary(scSide);

            if(entry.getValue() == EXON_INTRON)
            {
                // soft-clip must be intronic
                if(scSide == SE_START && readBoundary > regionBoundary)
                    continue;
                else if(scSide == SE_END && readBoundary < regionBoundary)
                    continue;

                boundaryDistance = abs(readBoundary - regionBoundary);

                if(boundaryDistance > MAX_INTRON_DISTANCE)
                    continue;
            }

            if(selectedTransData == null || !selectedTransData.IsCanonical && transData.IsCanonical)
            {
                selectedTransData = transData;
                spliceType = isSpliceAcceptor ? "acceptor" : "donor";
                exonBoundary = regionBoundary;
                exonBoundaryDistance = boundaryDistance;
                exonRank = transExonRef.ExonRank;
            }

            /*
            transcriptInfo.add(String.format("%s:%d:%s:%d:%d",
                    transData.TransName, transExonRef.ExonRank,
                        isSpliceAcceptor ? SPLICE_TYPE_ACCEPTOR : SPLICE_TYPE_DONOR, exonBoundary, exonBoundaryDistance));
             */
        }

        if(selectedTransData == null)
            return;

        String posKey = positionKey(read.orientation(), scSide, exonBoundary);

        // make note of any supplementary read (ie fusion candidate) and exclude any read matching one
        if(read.hasSuppAlignment())
        {
            mSupplementaryReadKeys.add(posKey);
            return;
        }

        int cohortFrequency = mCohortFrequency.getCohortFrequency(read.Chromosome, posKey);

        // for now, filter these
        //if(cohortFrequency > 0)
        //    return;

        final String geneId = selectedTransData.GeneId;
        GeneData geneData = mGenes.genes().stream()
                .filter(x -> x.GeneData.GeneId.equals(geneId))
                .map(x -> x.GeneData)
                .findFirst().orElse(null);

        // determine average base qual within the SC region
        int totalBaseQual = 0;
        int startIndex = scSide == SE_START ? 0 : read.Length - scLength;
        int endIndex = scSide == SE_START ? scLength : read.Length;
        for(int i = startIndex; i < endIndex; ++i)
        {
            totalBaseQual += record.getBaseQualities()[i];
        }

        double avgBaseQual = totalBaseQual / scLength;

        String scBases = scSide == SE_START ?
                Nucleotides.reverseStrandBases(read.ReadBases.substring(0, scLength)) : read.ReadBases.substring(read.ReadBases.length() - scLength);

        String mateCoords = String.format("%s:%d", read.mateChromosome(), read.mateStartPosition());

        UnmappedRead umRead = new UnmappedRead(
                read.Id, new ChrBaseRegion(read.Chromosome, read.PosStart, read.PosEnd), read.orientation(), scLength,
                scSide, avgBaseQual, geneData.GeneId, geneData.GeneName, selectedTransData.TransName, exonRank, exonBoundary,
                exonBoundaryDistance, spliceType, scBases, mateCoords, cohortFrequency, false);

        mCandidateReads.add(umRead);
    }

    public void writeUnmappedReads()
    {
        /*
        List<UnmappedRead> filteredReads = mCandidateReads.stream()
                .filter(x -> !mSupplementaryReadKeys.contains(x.positionKey()))
                .collect(Collectors.toList());

        writeReadData(mWriter, filteredReads);
        */

        mCandidateReads.forEach(x -> x.MatchesSupplementary = mSupplementaryReadKeys.contains(x.positionKey()));

        writeReadData(mWriter, mCandidateReads);
    }

    public static BufferedWriter createWriter(final IsofoxConfig config)
    {
        try
        {
            final String outputFileName = config.formOutputFile(UNMAPPED_READS_FILE_ID);

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write(UnmappedRead.header());
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to create unmapped reads writer: {}", e.toString());
            return null;
        }
    }

    private synchronized static void writeReadData(final BufferedWriter writer, final List<UnmappedRead> unmappedReads)
    {
        try
        {
            for(UnmappedRead umRead : unmappedReads)
            {
                writer.write(umRead.toCsv());
                writer.newLine();
            }
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write unmapped read data: {}", e.toString());
        }
    }

}