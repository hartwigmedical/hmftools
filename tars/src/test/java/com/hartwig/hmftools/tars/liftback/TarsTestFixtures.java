package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.liftback.rescue.RefSequenceSource;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

// Shared test fixtures for the liftback suite: the standard three-exon transcript contig, SAMRecord builders,
// a MockRefGenome-backed RefSequenceSource, and factories/builder for LiftedAlignment + LiftBackResult. Centralises
// what was previously copy-pasted into a dozen test classes (and keeps the giant LiftBackResult record off call sites).
public final class TarsTestFixtures
{
    public static final String GENOMIC_CHR = CHR_1;
    public static final String GENE_ID = "ENSG_TEST";
    public static final String GENE_NAME = "TESTG";
    public static final String TRANS_NAME = "ENST_TEST";
    public static final String TX_CONTIG = "ens" + GENE_ID + "_" + GENE_NAME + "_" + TRANS_NAME;

    private TarsTestFixtures() {}

    // exon spans on chr1: 100-199, 300-399, 500-549; introns 200-299, 400-499; contig length 250.
    public static ContigEntry threeExonContig()
    {
        return new ContigEntry(
                TX_CONTIG, 1, 250, GENE_ID, GENE_NAME, TRANS_NAME, GENOMIC_CHR, 1,
                List.of(new BaseRegion(100, 199), new BaseRegion(300, 399), new BaseRegion(500, 549)));
    }

    public static SAMRecord primaryRecord(final String contig, final int pos, final String cigar)
    {
        return primaryRecord("readX", contig, pos, cigar);
    }

    public static SAMRecord primaryRecord(final String readName, final String contig, final int pos, final String cigar)
    {
        final SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName(readName);
        record.setReferenceName(contig);
        record.setAlignmentStart(pos);
        record.setCigarString(cigar);
        record.setMappingQuality(60);
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(true);
        record.setProperPairFlag(true);
        return record;
    }

    public static SAMRecord secondMateRecord(final String contig, final int pos, final String cigar)
    {
        return secondMateRecord("readX", contig, pos, cigar);
    }

    public static SAMRecord secondMateRecord(final String readName, final String contig, final int pos, final String cigar)
    {
        final SAMRecord record = primaryRecord(readName, contig, pos, cigar);
        record.setFirstOfPairFlag(false);
        record.setSecondOfPairFlag(true);
        return record;
    }

    public static SAMRecord supplementaryRecord(final String contig, final int pos, final String cigar, final String saTag)
    {
        final SAMRecord record = supplementaryRecord("readX", contig, pos, cigar);
        record.setAttribute("SA", saTag);
        return record;
    }

    public static SAMRecord supplementaryRecord(final String readName, final String contig, final int pos, final String cigar)
    {
        final SAMRecord record = primaryRecord(readName, contig, pos, cigar);
        record.setSupplementaryAlignmentFlag(true);
        return record;
    }

    public static SAMRecord unpairedRecord(final String readName)
    {
        final SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName(readName);
        record.setReadPairedFlag(false);
        return record;
    }

    // RefSequenceSource over an in-memory chromosome (1-based, inclusive); empty out-of-range -> null per contract.
    public static RefSequenceSource refSource(final String chromosome, final String bases)
    {
        final MockRefGenome ref = new MockRefGenome(true);
        ref.RefGenomeMap.put(chromosome, bases);
        return (chrom, posStart, posEnd) ->
        {
            final byte[] result = ref.getBases(chrom, posStart, posEnd);
            return result.length == 0 ? null : result;
        };
    }

    public static byte[] repeatedBase(final int length, final char base)
    {
        final byte[] out = new byte[length];
        Arrays.fill(out, (byte) base);
        return out;
    }

    // ASCII bytes for a base string, e.g. read/ref windows in the tail-extend and collapse tests.
    public static byte[] bases(final String sequence)
    {
        return sequence.getBytes(StandardCharsets.US_ASCII);
    }

    // Mutable in-memory genome for the ref-dependent passes (tail-extend / collapse / canonicalize / rescue
    // ref-verify). Allocate a chromosome with with(), overwrite bases at 1-based coords with set(), then hand
    // asRefSource() to the engine. Out-of-range reads return null per the RefSequenceSource contract.
    public static final class TestGenome
    {
        private final Map<String,byte[]> mBases = new HashMap<>();

        public TestGenome with(final String chromosome, final int length, final char fill)
        {
            final byte[] seq = new byte[length];
            Arrays.fill(seq, (byte) fill);
            mBases.put(chromosome, seq);
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
            final byte[] seq = mBases.get(chromosome);
            for(int i = 0; i < sequence.length(); ++i)
                seq[oneBasedStart - 1 + i] = (byte) sequence.charAt(i);
            return this;
        }

        // overwrite a run of one base at a 1-based start (e.g. a divergent exon/intron stretch).
        public TestGenome set(final String chromosome, final int oneBasedStart, final int count, final char base)
        {
            final byte[] seq = mBases.get(chromosome);
            for(int i = 0; i < count; ++i)
                seq[oneBasedStart - 1 + i] = (byte) base;
            return this;
        }

        public RefSequenceSource asRefSource()
        {
            final MockRefGenome ref = new MockRefGenome(true);
            for(final Map.Entry<String,byte[]> entry : mBases.entrySet())
                ref.RefGenomeMap.put(entry.getKey(), new String(entry.getValue(), StandardCharsets.US_ASCII));
            return (chrom, posStart, posEnd) ->
            {
                final byte[] result = ref.getBases(chrom, posStart, posEnd);
                return result.length == 0 ? null : result;
            };
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
        return new LiftedAlignment(
                LiftedAlignment.AlignmentSource.SELF, chrom, pos, cigar,
                chrom, pos, cigar, 100, 0, null, null, null, false, true);
    }

    public static ResultBuilder resultBuilder()
    {
        return new ResultBuilder();
    }

    // Fluent builder over the 23-component LiftBackResult record with sensible defaults (a clean single-locus
    // ref-only primary). Tests override only the fields the scenario cares about.
    public static final class ResultBuilder
    {
        private LiftBackCategory mCategory = LiftBackCategory.REF_SINGLE;
        private LiftBackResult.Composition mComp = LiftBackResult.Composition.REF_ONLY;
        private LiftBackResult.RecordRole mRole = LiftBackResult.RecordRole.PRIMARY;
        private String mChrom = GENOMIC_CHR;
        private int mPos = 1000;
        private String mCigar = "50M";
        private boolean mNegativeStrand = false;
        private boolean mHasNCigar = false;
        private int mInputMapq = 60;
        private int mUpdatedMapq = 60;
        private int mNumXaAlts = 0;
        private int mNumRefAlts = 1;
        private int mNumTxAlts = 0;
        private int mNumLoci = 1;
        private int mNumDistinctCigarsAtPrimaryLocus = 1;
        private boolean mTxHasNCigar = false;
        private boolean mTxSoftClipAtBoundary = false;
        private boolean mRefSoftClipped = false;
        private boolean mRefFullMatch = true;
        private String mGeneIds = "";
        private String mNotes = "";
        private int mTranscriptStrand = 0;
        private List<LiftedAlignment> mAlignments = List.of();

        public ResultBuilder category(final LiftBackCategory v) { mCategory = v; return this; }
        public ResultBuilder comp(final LiftBackResult.Composition v) { mComp = v; return this; }
        public ResultBuilder role(final LiftBackResult.RecordRole v) { mRole = v; return this; }
        public ResultBuilder chrom(final String v) { mChrom = v; return this; }
        public ResultBuilder pos(final int v) { mPos = v; return this; }
        public ResultBuilder cigar(final String v) { mCigar = v; return this; }
        public ResultBuilder negativeStrand(final boolean v) { mNegativeStrand = v; return this; }
        public ResultBuilder hasNCigar(final boolean v) { mHasNCigar = v; return this; }
        public ResultBuilder inputMapq(final int v) { mInputMapq = v; return this; }
        public ResultBuilder updatedMapq(final int v) { mUpdatedMapq = v; return this; }
        public ResultBuilder numXaAlts(final int v) { mNumXaAlts = v; return this; }
        public ResultBuilder numRefAlts(final int v) { mNumRefAlts = v; return this; }
        public ResultBuilder numTxAlts(final int v) { mNumTxAlts = v; return this; }
        public ResultBuilder numLoci(final int v) { mNumLoci = v; return this; }
        public ResultBuilder numDistinctCigarsAtPrimaryLocus(final int v) { mNumDistinctCigarsAtPrimaryLocus = v; return this; }
        public ResultBuilder txHasNCigar(final boolean v) { mTxHasNCigar = v; return this; }
        public ResultBuilder txSoftClipAtBoundary(final boolean v) { mTxSoftClipAtBoundary = v; return this; }
        public ResultBuilder refSoftClipped(final boolean v) { mRefSoftClipped = v; return this; }
        public ResultBuilder refFullMatch(final boolean v) { mRefFullMatch = v; return this; }
        public ResultBuilder geneIds(final String v) { mGeneIds = v; return this; }
        public ResultBuilder notes(final String v) { mNotes = v; return this; }
        public ResultBuilder transcriptStrand(final int v) { mTranscriptStrand = v; return this; }
        public ResultBuilder alignments(final List<LiftedAlignment> v) { mAlignments = v; return this; }

        public LiftBackResult build()
        {
            return new LiftBackResult(
                    mCategory, mComp, mRole, mChrom, mPos, mCigar, mNegativeStrand, mHasNCigar,
                    mInputMapq, mUpdatedMapq, mNumXaAlts, mNumRefAlts, mNumTxAlts, mNumLoci,
                    mNumDistinctCigarsAtPrimaryLocus, mTxHasNCigar, mTxSoftClipAtBoundary, mRefSoftClipped,
                    mRefFullMatch, mGeneIds, mNotes, mTranscriptStrand, mAlignments);
        }
    }
}
