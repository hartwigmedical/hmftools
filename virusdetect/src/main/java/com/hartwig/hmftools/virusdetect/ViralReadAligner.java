package com.hartwig.hmftools.virusdetect;

import static java.lang.Math.ceil;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MISMATCHES_AND_DELETIONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.Streams.partitionStream;

import java.io.File;
import java.util.List;
import java.util.Objects;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.bwa.BwaMemAligner;
import com.hartwig.hmftools.common.bwa.BwaMemAlignParams;
import com.hartwig.hmftools.common.bwa.BwaMemAlignerConfig;
import com.hartwig.hmftools.common.bwa.IBwaMemAligner;
import com.hartwig.hmftools.common.codon.Nucleotides;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;

// Realigns candidate reads (a FASTA of read name + bases) to the viral reference with BWA-MEM in
// all-hits mode, writing the alignments to a BAM.
public class ViralReadAligner
{
    private final IBwaMemAligner mAligner;
    private final SAMFileHeader mHeader;
    private final double mMinAlignmentScoreFraction;
    private final int mChunkSize;

    private static final Logger LOGGER = LogManager.getLogger(ViralReadAligner.class);

    public static ViralReadAligner create(VirusConfig config, ViralReference reference)
    {
        if(!new File(config.viralBwaIndexImage()).exists())
        {
            throw new UserInputError("viral BWA index image not found: " + config.viralBwaIndexImage());
        }

        BwaMemAligner.initLibrary(config.bwaLibPath());
        IBwaMemAligner aligner = new BwaMemAligner(buildAlignerConfig(config));

        return new ViralReadAligner(
                aligner, buildHeader(reference.sequenceDictionary()), config.minAlignmentScoreFraction(), config.alignmentBatchSize());
    }

    ViralReadAligner(IBwaMemAligner aligner, SAMFileHeader header, double minAlignmentScoreFraction, int chunkSize)
    {
        mAligner = aligner;
        mHeader = header;
        mMinAlignmentScoreFraction = minAlignmentScoreFraction;
        mChunkSize = chunkSize;
    }

    public void align(String candidateFastaFile, String outputBamFile)
    {
        try(FastaSequenceFile fasta = new FastaSequenceFile(new File(candidateFastaFile), true);
                SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(mHeader, false, new File(outputBamFile)))
        {
            // Written a chunk at a time so peak heap is one chunk, not the whole sample.
            Stream<ReferenceSequence> reads = Stream.generate(fasta::nextSequence).takeWhile(Objects::nonNull);
            ChunkResult total = partitionStream(reads, mChunkSize)
                    .map(chunk -> alignChunk(chunk, writer))
                    .reduce(ChunkResult.EMPTY, ChunkResult::add);

            LOGGER.info("aligned {} candidate reads: {} with viral alignments, {} alignments written to {}",
                    total.totalReads(), total.alignedReads(), total.writtenAlignments(), outputBamFile);
        }
    }

    private ChunkResult alignChunk(List<ReferenceSequence> reads, SAMFileWriter writer)
    {
        List<byte[]> sequences = reads.stream().map(ReferenceSequence::getBases).toList();
        List<List<BwaMemAlignment>> alignments = mAligner.alignSequences(sequences);

        long alignedReads = 0;
        long writtenAlignments = 0;
        for(int i = 0; i < reads.size(); ++i)
        {
            List<SAMRecord> records = toRecords(reads.get(i).getName(), sequences.get(i), alignments.get(i));
            records.forEach(writer::addAlignment);
            if(!records.isEmpty())
            {
                ++alignedReads;
            }
            writtenAlignments += records.size();
        }
        return new ChunkResult(reads.size(), alignedReads, writtenAlignments);
    }

    // Alignments scoring below the read-length fraction, and reads with no hit, support no virus and are not written.
    private List<SAMRecord> toRecords(String readName, byte[] readBases, List<BwaMemAlignment> alignments)
    {
        int minScore = (int) ceil(mMinAlignmentScoreFraction * readBases.length);
        return alignments.stream()
                .filter(alignment -> alignment.getRefId() >= 0 && alignment.getAlignerScore() >= minScore)
                .map(alignment -> toRecord(readName, readBases, alignment))
                .toList();
    }

    // Built from the alignment alone: the FASTA carried only name and bases, so there are no base
    // qualities, no original tags, and no mate (single-end; pair mates differ only by name suffix).
    private SAMRecord toRecord(String readName, byte[] readBases, BwaMemAlignment alignment)
    {
        SAMRecord record = new SAMRecord(mHeader);
        record.setReadName(readName);
        record.setFlags(alignment.getSamFlag());
        record.setReferenceIndex(alignment.getRefId());
        record.setAlignmentStart(alignment.getRefStart() + 1);   // BWA reference positions are 0-based; SAM is 1-based
        record.setCigarString(alignment.getCigar());
        record.setMappingQuality(alignment.getMapQual());

        byte[] bases = record.getReadNegativeStrandFlag() ? Nucleotides.reverseComplementBases(readBases) : readBases;
        if(record.getCigar().getReadLength() == bases.length)
        {
            record.setReadBases(bases);
        }
        else
        {
            record.setReadBases(SAMRecord.NULL_SEQUENCE);   // hard-clipped alignment: sequence omitted, consistent with the CIGAR
        }

        record.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, alignment.getAlignerScore());
        record.setAttribute(NUM_MUTATONS_ATTRIBUTE, alignment.getNMismatches());
        record.setAttribute(MISMATCHES_AND_DELETIONS_ATTRIBUTE, alignment.getMDTag());
        return record;
    }

    private static BwaMemAlignerConfig buildAlignerConfig(VirusConfig config)
    {
        boolean allAlignments = true;
        return new BwaMemAlignerConfig(
                config.viralBwaIndexImage(), BwaMemAlignParams.DEFAULT, allAlignments, config.threads(), config.alignmentBatchSize());
    }

    private static SAMFileHeader buildHeader(SAMSequenceDictionary dictionary)
    {
        SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        return header;
    }

    private record ChunkResult(long totalReads, long alignedReads, long writtenAlignments)
    {
        static final ChunkResult EMPTY = new ChunkResult(0, 0, 0);

        ChunkResult add(ChunkResult other)
        {
            return new ChunkResult(
                    totalReads + other.totalReads, alignedReads + other.alignedReads, writtenAlignments + other.writtenAlignments);
        }
    }
}
