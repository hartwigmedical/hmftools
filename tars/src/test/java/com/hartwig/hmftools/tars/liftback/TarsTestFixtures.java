package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.tars.common.ContigEntry;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

// Shared test fixtures for the liftback suite: the standard three-exon transcript contig, SAMRecord builders,
// a MockRefGenome-backed reference genome, and factories/builder for LiftedAlignment + LiftedRecord.
public final class TarsTestFixtures
{
    public static final String GENOMIC_CHR = CHR_1;
    public static final String GENE_ID = "ENSG_TEST";
    public static final String GENE_NAME = "TESTG";
    public static final String TRANS_NAME = "ENST_TEST";
    public static final String TX_CONTIG = "ens" + GENE_ID + "_" + GENE_NAME + "_" + TRANS_NAME;

    private TarsTestFixtures() { }

    // exon spans on chr1: 100-199, 300-399, 500-549; introns 200-299, 400-499; contig length 250.
    public static ContigEntry threeExonContig()
    {
        return new ContigEntry(
                TX_CONTIG, 1, 250, GENE_ID, GENE_NAME, TRANS_NAME, GENOMIC_CHR, 1,
                List.of(new BaseRegion(100, 199), new BaseRegion(300, 399), new BaseRegion(500, 549)));
    }

    // an ExonRegionIndex from a single entry carrying the exon spans on the given chromosome (bare "1" or "chr1" both
    // normalize to the V38 chr1 key the lift emits).
    public static ExonRegionIndex exonRegionIndex(final String chromosome, final List<int[]> exonSpans)
    {
        List<BaseRegion> spans = new ArrayList<>();
        for(int[] span : exonSpans)
        {
            spans.add(new BaseRegion(span[0], span[1]));
        }
        ContigEntry entry = ContigEntry.annotationOnly(GENE_ID, GENE_NAME, TRANS_NAME, V38.versionedChromosome(chromosome), 1, spans);
        return ExonRegionIndex.fromContigEntries(List.of(entry));
    }

    public static SAMRecord primaryRecord(final String contig, final int pos, final String cigar)
    {
        return primaryRecord("readX", contig, pos, cigar);
    }

    public static SAMRecord primaryRecord(final String readName, final String contig, final int pos, final String cigar)
    {
        SAMRecord record = mappedRecord(readName, contig, pos, cigar);
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(true);
        record.setProperPairFlag(true);
        return record;
    }

    // Mapped record with no pairing flags set - the single-end (Ultima) shape, and the base every other builder layers on.
    public static SAMRecord mappedRecord(final String readName, final String contig, final int pos, final String cigar)
    {
        SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName(readName);
        record.setReferenceName(contig);
        record.setAlignmentStart(pos);
        record.setCigarString(cigar);
        record.setMappingQuality(60);
        return record;
    }

    // Mapped record on either side of a pair, with an explicit strand - for mate-field patching.
    public static SAMRecord pairedRecord(
            final String readName, final boolean firstOfPair, final String contig, final int pos, final String cigar,
            final boolean negativeStrand)
    {
        SAMRecord record = mappedRecord(readName, contig, pos, cigar);
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(firstOfPair);
        record.setSecondOfPairFlag(!firstOfPair);
        record.setReadNegativeStrandFlag(negativeStrand);
        return record;
    }

    public static SAMRecord unmappedRecord(final String readName)
    {
        SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName(readName);
        record.setReadUnmappedFlag(true);
        return record;
    }

    public static SAMRecord pairedUnmappedRecord(final String readName, final boolean firstOfPair)
    {
        SAMRecord record = unmappedRecord(readName);
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(firstOfPair);
        record.setSecondOfPairFlag(!firstOfPair);
        return record;
    }

    public static SAMRecord secondMateRecord(final String contig, final int pos, final String cigar)
    {
        return secondMateRecord("readX", contig, pos, cigar);
    }

    public static SAMRecord secondMateRecord(final String readName, final String contig, final int pos, final String cigar)
    {
        SAMRecord record = primaryRecord(readName, contig, pos, cigar);
        record.setFirstOfPairFlag(false);
        record.setSecondOfPairFlag(true);
        return record;
    }

    public static SAMRecord supplementaryRecord(final String contig, final int pos, final String cigar, final String saTag)
    {
        SAMRecord record = supplementaryRecord("readX", contig, pos, cigar);
        record.setAttribute("SA", saTag);
        return record;
    }

    public static SAMRecord supplementaryRecord(final String readName, final String contig, final int pos, final String cigar)
    {
        SAMRecord record = primaryRecord(readName, contig, pos, cigar);
        record.setSupplementaryAlignmentFlag(true);
        return record;
    }

    public static SAMRecord unpairedRecord(final String readName)
    {
        SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName(readName);
        record.setReadPairedFlag(false);
        return record;
    }

    // mapped single-end (Ultima) primary: like primaryRecord but with no pairing flags set at all.
    public static SAMRecord unpairedPrimaryRecord(final String contig, final int pos, final String cigar)
    {
        return unpairedPrimaryRecord("readX", contig, pos, cigar);
    }

    public static SAMRecord unpairedPrimaryRecord(final String readName, final String contig, final int pos, final String cigar)
    {
        SAMRecord record = mappedRecord(readName, contig, pos, cigar);
        record.setReadPairedFlag(false);
        return record;
    }

    public static SAMRecord unpairedSupplementaryRecord(final String contig, final int pos, final String cigar, final String saTag)
    {
        SAMRecord record = unpairedPrimaryRecord(contig, pos, cigar);
        record.setSupplementaryAlignmentFlag(true);
        record.setAttribute("SA", saTag);
        return record;
    }

    // reference genome over a single in-memory chromosome
    public static RefGenomeInterface refGenome(final String chromosome, final String bases)
    {
        return new TestGenome().with(chromosome, bases(bases)).asRefGenome();
    }

    public static byte[] repeatedBase(final int length, final char base)
    {
        byte[] out = new byte[length];
        Arrays.fill(out, (byte) base);
        return out;
    }

    public static byte[] bases(final String sequence)
    {
        return sequence.getBytes(StandardCharsets.US_ASCII);
    }

    // Mutable in-memory genome for the ref-dependent passes. Allocate a chromosome with with(), overwrite bases at
    // 1-based coords with set(), then hand asRefGenome() to the engine.
    public static final class TestGenome
    {
        private final Map<String, byte[]> mBases = new HashMap<>();

        public TestGenome with(final String chromosome, final int length, final char fill)
        {
            mBases.put(chromosome, repeatedBase(length, fill));
            return this;
        }

        public TestGenome with(final String chromosome, final byte[] bases)
        {
            mBases.put(chromosome, bases);
            return this;
        }

        // overwrite bases at a 1-based start (e.g. to match a read for NM=0, or seed a splice motif).
        public TestGenome set(final String chromosome, final int oneBasedStart, final String sequence)
        {
            byte[] seq = mBases.get(chromosome);
            for(int i = 0; i < sequence.length(); ++i)
            {
                seq[oneBasedStart - 1 + i] = (byte) sequence.charAt(i);
            }
            return this;
        }

        // overwrite a run of one base at a 1-based start (e.g. a divergent exon/intron stretch).
        public TestGenome set(final String chromosome, final int oneBasedStart, final int count, final char base)
        {
            byte[] seq = mBases.get(chromosome);
            for(int i = 0; i < count; ++i)
            {
                seq[oneBasedStart - 1 + i] = (byte) base;
            }
            return this;
        }

        public RefGenomeInterface asRefGenome()
        {
            MockRefGenome refGenome = new MockRefGenome(true);
            for(Map.Entry<String, byte[]> entry : mBases.entrySet())
            {
                refGenome.RefGenomeMap.put(entry.getKey(), new String(entry.getValue(), StandardCharsets.US_ASCII));
            }
            return refGenome;
        }
    }

    // assert a lifted record's genomic placement (chrom + start + cigar) in one line.
    public static void assertLifted(final SAMRecord record, final String chrom, final int pos, final String cigar)
    {
        assertEquals(chrom, record.getReferenceName());
        assertEquals(pos, record.getAlignmentStart());
        assertEquals(cigar, record.getCigarString());
    }

    public static LiftedAlignment selfAlignment(final String chrom, final int pos, final String cigar)
    {
        return new LiftedAlignment(chrom, pos, cigar, 0, false, false, true, 0);
    }

    // A lifted record holding a single placement - the shape a mate is read as.
    public static LiftedRecord liftedRecordAt(
            final String chromosome, final int pos, final String cigar, final boolean negativeStrand)
    {
        LiftedAlignment alignment = new LiftedAlignment(chromosome, pos, cigar, 0, false, false, !negativeStrand, 0);
        return recordBuilder().alignments(List.of(alignment)).build();
    }

    public static RecordBuilder recordBuilder()
    {
        return new RecordBuilder();
    }

    // Fluent builder over LiftedRecord, defaulting to a clean single-locus primary. The placement comes from the
    // supplied alignments (primaryIndex, default 0), so tests set only the decision fields they care about.
    public static final class RecordBuilder
    {
        private int mUpdatedMapQuality = 60;
        private int mNumLoci = 1;
        private String mNotes = "";
        private int mPrimaryIndex = 0;
        private List<LiftedAlignment> mAlignments = List.of();

        public RecordBuilder updatedMapQuality(final int value)
        {
            mUpdatedMapQuality = value;
            return this;
        }

        public RecordBuilder numLoci(final int value)
        {
            mNumLoci = value;
            return this;
        }

        public RecordBuilder notes(final String value)
        {
            mNotes = value;
            return this;
        }

        public RecordBuilder primaryIndex(final int value)
        {
            mPrimaryIndex = value;
            return this;
        }

        public RecordBuilder alignments(final List<LiftedAlignment> value)
        {
            mAlignments = value;
            return this;
        }

        public LiftedRecord build()
        {
            return new LiftedRecord(mUpdatedMapQuality, mNumLoci, mNotes, mPrimaryIndex, mAlignments);
        }
    }
}
