package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.common.TarsConstants.MAX_IMPLIED_INTRON_LENGTH;
import static com.hartwig.hmftools.tars.common.TarsConstants.MAX_SUPP_MERGES;
import static com.hartwig.hmftools.tars.common.TarsConstants.MAX_SUPP_READ_OVERLAP;
import static com.hartwig.hmftools.tars.common.TarsConstants.MIN_IMPLIED_INTRON_LENGTH;

import static org.junit.Assert.assertEquals;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.liftback.features.SupplementaryConfig;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

// Shared fixtures for the liftback suite: the three-exon transcript contig, SAMRecord builders, an in-memory reference
// genome and factories for LiftedAlignment + LiftedRecord.
public final class TarsTestFixtures
{
    public static final String GENOMIC_CHR = CHR_1;
    public static final String GENE_ID = "ENSG_TEST";
    public static final String GENE_NAME = "TESTG";
    public static final String TRANS_NAME = "ENST_TEST";
    public static final String TX_CONTIG = "ens" + GENE_ID + "_" + GENE_NAME + "_" + TRANS_NAME;

    private TarsTestFixtures() { }

    public static SupplementaryConfig supplementaryConfig()
    {
        return new SupplementaryConfig(
                MIN_IMPLIED_INTRON_LENGTH, MAX_IMPLIED_INTRON_LENGTH,
                MAX_SUPP_MERGES, false, MAX_SUPP_READ_OVERLAP);
    }

    // exon spans on chr1: 100-199, 300-399, 500-549; introns 200-299, 400-499; contig length 250.
    public static ContigEntry threeExonContig()
    {
        return new ContigEntry(
                TX_CONTIG, 1, 250, GENE_ID, GENE_NAME, TRANS_NAME, GENOMIC_CHR, 1,
                List.of(new BaseRegion(100, 199), new BaseRegion(300, 399), new BaseRegion(500, 549)));
    }

    // bare "1" and "chr1" both normalize to the V38 chr1 key the lift emits.
    public static EnsemblAnnotationIndex exonRegionIndex(final String chromosome, final List<int[]> exonSpans)
    {
        List<BaseRegion> spans = new ArrayList<>();
        for(int[] span : exonSpans)
        {
            spans.add(new BaseRegion(span[0], span[1]));
        }
        ContigEntry entry = ContigEntry.annotationOnly(GENE_ID, GENE_NAME, TRANS_NAME, V38.versionedChromosome(chromosome), 1, spans);
        return EnsemblAnnotationIndex.fromContigEntries(List.of(entry));
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

    // Mapped record on either side of a pair with an explicit strand - for mate-field patching.
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

    // Mapped single-end (Ultima) primary: like primaryRecord but with no pairing flags set at all.
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

    // Reference genome over a single in-memory chromosome.
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

    // Mutable in-memory genome for the ref-dependent passes; set() coordinates are 1-based.
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

        // Overwrite bases to match a read for NM=0, or to seed a splice motif.
        public TestGenome set(final String chromosome, final int oneBasedStart, final String sequence)
        {
            byte[] seq = mBases.get(chromosome);
            for(int i = 0; i < sequence.length(); ++i)
            {
                seq[oneBasedStart - 1 + i] = (byte) sequence.charAt(i);
            }
            return this;
        }

        // Overwrite a run of one base, e.g. a divergent exon/intron stretch.
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

    // Assert a lifted record's genomic placement: chrom, start and cigar.
    public static void assertLifted(final SAMRecord record, final String chrom, final int pos, final String cigar)
    {
        assertEquals(chrom, record.getReferenceName());
        assertEquals(pos, record.getAlignmentStart());
        assertEquals(cigar, record.getCigarString());
    }

    public static LiftedAlignment selfAlignment(final String chrom, final int pos, final String cigar)
    {
        return new LiftedAlignment(chrom, pos, cigar, 0, false, true, 0);
    }

    // Single-placement lifted record - the shape a mate is read as.
    public static LiftedRecord liftedRecordAt(
            final String chromosome, final int pos, final String cigar, final boolean negativeStrand)
    {
        LiftedAlignment alignment = new LiftedAlignment(chromosome, pos, cigar, 0, false, !negativeStrand, 0);
        return recordBuilder().alignments(List.of(alignment)).build();
    }

    public static RecordBuilder recordBuilder()
    {
        return new RecordBuilder();
    }

    // Builds a LiftedRecord defaulting to a clean single-locus primary; placement comes from the supplied alignments at primaryIndex.
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
