package com.hartwig.hmftools.isofox.unmapped;

import static java.lang.Math.abs;
import static java.lang.Math.scalb;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxFunction.UNMAPPED_READS;
import static com.hartwig.hmftools.isofox.common.ReadRecord.findOverlappingRegions;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_BOUNDARY;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_INTRON;
import static com.hartwig.hmftools.isofox.unmapped.UnmappedRead.UMR_NO_MATE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.bam.ClippedSide;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.FragmentTracker;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.common.TransExonRef;

public class UmrFinder
{
    private final IsofoxConfig mConfig;
    private final boolean mEnabled;
    private GeneCollection mGenes;

    private final Set<String> mChimericReadKeys;
    private final List<UnmappedRead> mCandidateReads;

    private final BufferedWriter mWriter;

    public static final String UNMAPPED_READS_FILE_ID = "unmapped_reads.csv";

    private static final int MIN_SOFT_CLIP_LENGTH = 11;
    private static final int MAX_INTRON_DISTANCE = 5;

    public UmrFinder(final IsofoxConfig config, final BufferedWriter writer)
    {
        mEnabled = config.runFunction(UNMAPPED_READS);
        mConfig = config;

        // mCandidateReads = Lists.newArrayList();
        mChimericReadKeys = Sets.newHashSet();
        mCandidateReads = Lists.newArrayList();
        mGenes = null;

        mWriter = writer;
    }

    public boolean enabled() { return mEnabled; }

    public Set<String> getSupplementaryReadKeys() { return mChimericReadKeys; }
    public List<UnmappedRead> getCandidateReads() { return mCandidateReads; }

    public void setGeneData(final GeneCollection genes)
    {
        mGenes = genes;
        mCandidateReads.clear();
        mChimericReadKeys.clear();
    }

    public void processReads(final ReadRecord read1, final ReadRecord read2, boolean isChimeric)
    {
        if(!mEnabled)
            return;

        // check chimeric reads for soft-clips on or near an exon boundary to exclude these from the candidate unmapped collection
        UnmappedRead umRead1 = evaluateRead(read1);
        UnmappedRead umRead2 = evaluateRead(read2);

        if(umRead1 == null && umRead2 == null)
            return;

        // make note of any supplementary read (ie fusion candidate) and exclude any read matching one
        if(isChimeric)
        {
            if(umRead1 != null)
                mChimericReadKeys.add(umRead1.positionKey());

            if(umRead2 != null)
                mChimericReadKeys.add(umRead2.positionKey());

            return;
        }

        // exclude the unmapped read(s) if:
        // a) they have soft-clips on both sides
        // b) one read is past the unmapped soft-clip side

        if(umRead1 != null && umRead2 != null && umRead1.matches(umRead2))
        {
            // both reads support the same junction
        }
        else
        {
            if(umRead1 != null && !isValidReadPair(umRead1, read2))
                return;

            if(umRead2 != null && !isValidReadPair(umRead2, read1))
                return;
        }

        if(umRead1 != null)
            mCandidateReads.add(umRead1);

        if(umRead2 != null)
            mCandidateReads.add(umRead2);
    }

    private boolean isValidReadPair(final UnmappedRead umRead, final ReadRecord otherRead)
    {
        int otherScSidePosition = otherRead.getCoordsBoundary(umRead.ScSide);

        if(umRead.ScSide == SE_START && otherScSidePosition < umRead.ReadRegion.start())
            return false;

        if(umRead.ScSide == SE_END && otherScSidePosition > umRead.ReadRegion.end())
            return false;

        // ignore reads where they overlap and their fragment size is less than the read length
        if(otherRead.likelyAdaperSoftClipping())
            return false;

        if(otherRead.isSoftClipped(switchIndex(umRead.ScSide)))
        {
            int scLength = umRead.ScSide == SE_START ?
                    otherRead.Cigar.getLastCigarElement().getLength() : otherRead.Cigar.getFirstCigarElement().getLength();

            return scLength <= 5; // all small SCs only
        }

        return true;
    }

    public void processUnpairedReads(final FragmentTracker fragmentTracker)
    {
        for(Object object : fragmentTracker.getValues())
        {
            final ReadRecord read = (ReadRecord)object;

            UnmappedRead umRead = evaluateRead(read);

            if(umRead != null)
            {
                if(read.isChimeric())
                    mChimericReadKeys.add(umRead.positionKey());
                else
                    mCandidateReads.add(umRead);
            }
        }
    }

    public UnmappedRead evaluateRead(final ReadRecord read)
    {
        // find reads with sufficient soft-clipping and overlapping a splice junction
        if(!read.containsSoftClipping())
            return null;

        ClippedSide clippedSide = ClippedSide.fromCigar(read.Cigar, false);

        if(clippedSide == null)
            return null;

        int scSide = clippedSide.Side;
        int scLength = clippedSide.Length;

        if(scLength < MIN_SOFT_CLIP_LENGTH)
            return null;

        final List<RegionReadData> overlappingRegions = findOverlappingRegions(mGenes.getExonRegions(), read);

        if(overlappingRegions.isEmpty())
            return null;

        read.processOverlappingRegions(overlappingRegions);

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

            // ignore start and end of transcript
            if(isSpliceAcceptor && transExonRef.ExonRank == 1)
                continue;

            if(!isSpliceAcceptor && transExonRef.ExonRank == transData.exons().size())
                continue;

            int regionBoundary = scSide == SE_START ? region.start() : region.end();
            int readBoundary = read.getCoordsBoundary(scSide);

            // soft-clip must be intronic if the read spans into the intron
            if(scSide == SE_START && readBoundary > regionBoundary)
                continue;
            else if(scSide == SE_END && readBoundary < regionBoundary)
                continue;

            int boundaryDistance = abs(readBoundary - regionBoundary);

            if(boundaryDistance > MAX_INTRON_DISTANCE)
                continue;

            if(selectedTransData == null || !selectedTransData.IsCanonical && transData.IsCanonical)
            {
                selectedTransData = transData;
                spliceType = isSpliceAcceptor ? "acceptor" : "donor";
                exonBoundary = regionBoundary;
                exonBoundaryDistance = boundaryDistance;
                exonRank = transExonRef.ExonRank;
            }
        }

        if(selectedTransData == null)
            return null;

        final String geneId = selectedTransData.GeneId;
        GeneData geneData = mGenes.genes().stream()
                .filter(x -> x.GeneData.GeneId.equals(geneId))
                .map(x -> x.GeneData)
                .findFirst().orElse(null);

        // determine average base qual within the SC region
        double avgBaseQual = 0;

        if(read.baseQualities() != null)
        {
            int totalBaseQual = 0;
            int startIndex = scSide == SE_START ? 0 : read.Length - scLength;
            int endIndex = scSide == SE_START ? scLength : read.Length;
            for(int i = startIndex; i < endIndex; ++i)
            {
                totalBaseQual += read.baseQualities()[i];
            }

            avgBaseQual = totalBaseQual / scLength;
        }

        // take all SC bases but also any from the exon junction
        int scBaseLength = scLength + exonBoundaryDistance;

        String scBases = scSide == SE_START ?
                Nucleotides.reverseComplementBases(read.ReadBases.substring(0, scBaseLength)) : read.ReadBases.substring(read.ReadBases.length() - scBaseLength);

        byte[] scBaseQuals = new byte[scBaseLength];

        if(read.baseQualities() != null)
        {
            if(scSide == SE_START)
            {
                for(int i = 0; i < scBaseLength; ++i)
                {
                    int index = scBaseLength - i - 1;
                    scBaseQuals[i] = read.baseQualities()[index];
                }
            }
            else
            {
                int bqLength = read.baseQualities().length;
                int offset = bqLength - scBaseLength;
                for(int i = 0; i < scBaseLength; ++i)
                {
                    int index = i + offset;
                    scBaseQuals[i] = read.baseQualities()[index];
                }
            }
        }

        if(scBases.contains("N"))
            return null;

        String mateCoords = read.isMateUnmapped() ?
                UMR_NO_MATE : String.format("%s:%d", read.mateChromosome(), read.mateStartPosition());

        UnmappedRead umRead = new UnmappedRead(
                read.Id, new ChrBaseRegion(read.Chromosome, read.PosStart, read.PosEnd), scLength,
                scSide, avgBaseQual, geneData.GeneId, geneData.GeneName, selectedTransData.TransName, exonRank, exonBoundary,
                exonBoundaryDistance, spliceType, scBases, scBaseQuals, mateCoords, false);

        return umRead;
    }

    public void markChimericMatches()
    {
        mCandidateReads.forEach(x -> x.MatchesChimeric = mChimericReadKeys.contains(x.positionKey()));
    }

    public void writeUnmappedReads()
    {
        markChimericMatches();

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
        if(writer == null)
            return;

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