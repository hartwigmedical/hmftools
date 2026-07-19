package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createGeneDataCache;

import static org.junit.Assert.assertEquals;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.liftback.supplementary.RefSequenceSource;

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

    private TarsTestFixtures() { }

    // exon spans on chr1: 100-199, 300-399, 500-549; introns 200-299, 400-499; contig length 250.
    public static ContigEntry threeExonContig()
    {
        return new ContigEntry(
                TX_CONTIG, 1, 250, GENE_ID, GENE_NAME, TRANS_NAME, GENOMIC_CHR, 1,
                List.of(new BaseRegion(100, 199), new BaseRegion(300, 399), new BaseRegion(500, 549)));
    }

    // builds an ExonRegionIndex from a single gene on the given ensembl chromosome (bare, e.g. "1") whose one
    // transcript carries the exon spans; fromCache merges the spans and keys them in V38 form (chr1). Shared by the
    // resolver and exon-index suites.
    public static ExonRegionIndex exonRegionIndex(final String chromosome, final List<int[]> exonSpans)
    {
        EnsemblDataCache cache = createGeneDataCache();
        addGeneData(cache, chromosome, List.of(createEnsemblGeneData(GENE_ID, GENE_NAME, chromosome, 1, 1, 100_000)));

        TranscriptData transcript = new TranscriptData(
                1, TRANS_NAME, GENE_ID, true, (byte) 1, 1, 100_000, null, null, "protein_coding", "");
        List<ExonData> exons = new ArrayList<>();
        int rank = 1;
        for(int[] span : exonSpans)
        {
            exons.add(new ExonData(1, span[0], span[1], rank++, -1, -1));
        }
        transcript.setExons(exons);
        addTransExonData(cache, GENE_ID, List.of(transcript));

        return ExonRegionIndex.fromCache(cache, V38);
    }

    public static SAMRecord primaryRecord(final String contig, final int pos, final String cigar)
    {
        return primaryRecord("readX", contig, pos, cigar);
    }

    public static SAMRecord primaryRecord(final String readName, final String contig, final int pos, final String cigar)
    {
        SAMRecord record = new SAMRecord(new SAMFileHeader());
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
        SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName("readX");
        record.setReferenceName(contig);
        record.setAlignmentStart(pos);
        record.setCigarString(cigar);
        record.setMappingQuality(60);
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

    // RefSequenceSource over an in-memory chromosome (1-based, inclusive); empty out-of-range -> null per contract.
    public static RefSequenceSource refSource(final String chromosome, final String bases)
    {
        MockRefGenome ref = new MockRefGenome(true);
        ref.RefGenomeMap.put(chromosome, bases);
        return (chrom, posStart, posEnd) ->
        {
            byte[] result = ref.getBases(chrom, posStart, posEnd);
            return result.length == 0 ? null : result;
        };
    }

    public static byte[] repeatedBase(final int length, final char base)
    {
        byte[] out = new byte[length];
        Arrays.fill(out, (byte) base);
        return out;
    }

    // ASCII bytes for a base string, e.g. read/ref windows in the tail-extend and collapse tests.
    public static byte[] bases(final String sequence)
    {
        return sequence.getBytes(StandardCharsets.US_ASCII);
    }

    // Mutable in-memory genome for the ref-dependent passes (tail-extend / collapse / canonicalize / supplementary
    // resolve ref-verify). Allocate a chromosome with with(), overwrite bases at 1-based coords with set(), then hand
    // asRefSource() to the engine. Out-of-range reads return null per the RefSequenceSource contract.
    public static final class TestGenome
    {
        private final Map<String, byte[]> mBases = new HashMap<>();

        public TestGenome with(final String chromosome, final int length, final char fill)
        {
            byte[] seq = new byte[length];
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

        public RefSequenceSource asRefSource()
        {
            MockRefGenome ref = new MockRefGenome(true);
            for(Map.Entry<String, byte[]> entry : mBases.entrySet())
            {
                ref.RefGenomeMap.put(entry.getKey(), new String(entry.getValue(), StandardCharsets.US_ASCII));
            }
            return (chrom, posStart, posEnd) ->
            {
                byte[] result = ref.getBases(chrom, posStart, posEnd);
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
        private RecordState mRecordState = RecordState.RESOLVED;
        private DecidingFeature mDecidingFeature = DecidingFeature.SOLE_REF;
        private boolean mSwapped = false;
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

        public ResultBuilder recordState(final RecordState value)
        {
            mRecordState = value;
            return this;
        }

        public ResultBuilder decidingFeature(final DecidingFeature value)
        {
            mDecidingFeature = value;
            return this;
        }

        public ResultBuilder swapped(final boolean value)
        {
            mSwapped = value;
            return this;
        }

        public ResultBuilder comp(final LiftBackResult.Composition value)
        {
            mComp = value;
            return this;
        }

        public ResultBuilder role(final LiftBackResult.RecordRole value)
        {
            mRole = value;
            return this;
        }

        public ResultBuilder chrom(final String value)
        {
            mChrom = value;
            return this;
        }

        public ResultBuilder pos(final int value)
        {
            mPos = value;
            return this;
        }

        public ResultBuilder cigar(final String value)
        {
            mCigar = value;
            return this;
        }

        public ResultBuilder negativeStrand(final boolean value)
        {
            mNegativeStrand = value;
            return this;
        }

        public ResultBuilder hasNCigar(final boolean value)
        {
            mHasNCigar = value;
            return this;
        }

        public ResultBuilder inputMapq(final int value)
        {
            mInputMapq = value;
            return this;
        }

        public ResultBuilder updatedMapq(final int value)
        {
            mUpdatedMapq = value;
            return this;
        }

        public ResultBuilder numXaAlts(final int value)
        {
            mNumXaAlts = value;
            return this;
        }

        public ResultBuilder numRefAlts(final int value)
        {
            mNumRefAlts = value;
            return this;
        }

        public ResultBuilder numTxAlts(final int value)
        {
            mNumTxAlts = value;
            return this;
        }

        public ResultBuilder numLoci(final int value)
        {
            mNumLoci = value;
            return this;
        }

        public ResultBuilder numDistinctCigarsAtPrimaryLocus(final int value)
        {
            mNumDistinctCigarsAtPrimaryLocus = value;
            return this;
        }

        public ResultBuilder txHasNCigar(final boolean value)
        {
            mTxHasNCigar = value;
            return this;
        }

        public ResultBuilder txSoftClipAtBoundary(final boolean value)
        {
            mTxSoftClipAtBoundary = value;
            return this;
        }

        public ResultBuilder refSoftClipped(final boolean value)
        {
            mRefSoftClipped = value;
            return this;
        }

        public ResultBuilder refFullMatch(final boolean value)
        {
            mRefFullMatch = value;
            return this;
        }

        public ResultBuilder geneIds(final String value)
        {
            mGeneIds = value;
            return this;
        }

        public ResultBuilder notes(final String value)
        {
            mNotes = value;
            return this;
        }

        public ResultBuilder transcriptStrand(final int value)
        {
            mTranscriptStrand = value;
            return this;
        }

        public ResultBuilder alignments(final List<LiftedAlignment> value)
        {
            mAlignments = value;
            return this;
        }

        public LiftBackResult build()
        {
            return new LiftBackResult(
                    mRecordState, mDecidingFeature, mSwapped, mComp, mRole, mChrom, mPos, mCigar, mNegativeStrand, mHasNCigar,
                    mInputMapq, mUpdatedMapq, mNumXaAlts, mNumRefAlts, mNumTxAlts, mNumLoci,
                    mNumDistinctCigarsAtPrimaryLocus, mTxHasNCigar, mTxSoftClipAtBoundary, mRefSoftClipped,
                    mRefFullMatch, mGeneIds, mNotes, mTranscriptStrand, mAlignments);
        }
    }
}
