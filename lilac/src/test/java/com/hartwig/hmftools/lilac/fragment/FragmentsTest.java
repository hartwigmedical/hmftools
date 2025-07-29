package com.hartwig.hmftools.lilac.fragment;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.NANOS_IN_MILLISECOND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_BASE_QUAL;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.calcAminoAcidIndices;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.mergeFragments;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_A;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_B;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_C;
import static com.hartwig.hmftools.lilac.misc.LilacTestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.lilac.misc.LilacTestUtils.createReadRecord;
import static com.hartwig.hmftools.lilac.read.Read.createRead;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.lilac.read.Read;

import org.junit.Ignore;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class FragmentsTest
{
    private static final String TEST_READ_BASES = MockRefGenome.generateRandomBases(100);

    @Test
    public void testIndices()
    {
        assertRange(0, -1, calcAminoAcidIndices(0, 0));
        assertRange(0, -1, calcAminoAcidIndices(0, 1));
        assertRange(0, 0, calcAminoAcidIndices(0, 2));
        assertRange(1, 0, calcAminoAcidIndices(1, 2));
        assertRange(1, 0, calcAminoAcidIndices(2, 2));

        assertRange(1, 1, calcAminoAcidIndices(3, 5));
        assertRange(1, 1, calcAminoAcidIndices(3, 6));
        assertRange(1, 1, calcAminoAcidIndices(3, 7));
        assertRange(1, 2, calcAminoAcidIndices(3, 8));
        assertRange(1, 2, calcAminoAcidIndices(3, 9));
        assertRange(1, 2, calcAminoAcidIndices(3, 10));
    }

    @Test
    public void testReadPositions()
    {
        BaseRegion codingRegion = new BaseRegion(100, 200);

        // test the low-qual bases are trimming up until the lowest score point only
        SAMRecord record = SamRecordTestUtils.createSamRecord(
                TEST_READ_ID, CHR_1, 120, TEST_READ_BASES.substring(0, 50), "10S30M10S", CHR_1, 300,
                false, false, null);

        Read read = createRead(codingRegion, record, false, false);
        assertNotNull(read);

        assertEquals(120, read.PositionStart);
        assertEquals(149, read.PositionEnd);
        assertEquals(0, read.SoftClippedStart);
        assertEquals(0, read.SoftClippedEnd);
        assertEquals(10, read.ReadIndexStart);
        assertEquals(39, read.ReadIndexEnd);

        read = createRead(codingRegion, record, true, false);
        assertNotNull(read);

        assertEquals(110, read.PositionStart);
        assertEquals(159, read.PositionEnd);
        assertEquals(10, read.SoftClippedStart);
        assertEquals(10, read.SoftClippedEnd);
        assertEquals(0, read.ReadIndexStart);
        assertEquals(49, read.ReadIndexEnd);

        // now capped by the coding region
        codingRegion = new BaseRegion(115, 155);

        read = createRead(codingRegion, record, true, false);
        assertNotNull(read);

        assertEquals(115, read.PositionStart);
        assertEquals(155, read.PositionEnd);
        assertEquals(5, read.SoftClippedStart);
        assertEquals(6, read.SoftClippedEnd);
        assertEquals(5, read.ReadIndexStart);
        assertEquals(45, read.ReadIndexEnd);

        // within the aligned positions
        codingRegion = new BaseRegion(125, 135);

        read = createRead(codingRegion, record, true, false);
        assertNotNull(read);

        assertEquals(125, read.PositionStart);
        assertEquals(135, read.PositionEnd);
        assertEquals(0, read.SoftClippedStart);
        assertEquals(0, read.SoftClippedEnd);
        assertEquals(15, read.ReadIndexStart);
        assertEquals(25, read.ReadIndexEnd);

    }

    @Test
    public void testReadLowQualBaseTrimming()
    {
        BaseRegion codingRegion = new BaseRegion(100, 1000);

        // test the low-qual bases are trimming up until the lowest score point only
        SAMRecord record = SamRecordTestUtils.createSamRecord(
                TEST_READ_ID, CHR_1, 100, TEST_READ_BASES.substring(0, 50), "50M", CHR_1, 300,
                false, false, null);

        byte lowBaseQual = DEFAULT_MIN_BASE_QUAL - 1;
        setBaseQualities(record, 40, 49, lowBaseQual);

        Read read = createRead(codingRegion, record, true, false);
        assertNotNull(read);

        assertEquals(10, read.trimmedBases());
        assertEquals(100, read.PositionStart);
        assertEquals(139, read.PositionEnd);
        assertEquals(0, read.ReadIndexStart);
        assertEquals(39, read.ReadIndexEnd);

        // with a second low point but less than the first
        setBaseQualities(record, 28, 31, lowBaseQual);

        read = createRead(codingRegion, record, false, false);
        assertNotNull(read);

        assertEquals(10, read.trimmedBases());
        assertEquals(100, read.PositionStart);
        assertEquals(139, read.PositionEnd);
        assertEquals(0, read.ReadIndexStart);
        assertEquals(39, read.ReadIndexEnd);

        // test that trimmed bases factor in the coding region - first low qual ending mid way between read end and coding end
        codingRegion = new BaseRegion(100, 130);

        read = createRead(codingRegion, record, false, false);
        assertNotNull(read);

        assertEquals(10, read.trimmedBases());
        assertEquals(100, read.PositionStart);
        assertEquals(130, read.PositionEnd);
        assertEquals(0, read.ReadIndexStart);
        assertEquals(30, read.ReadIndexEnd);

        // now low qual extending inside the coding region
        setBaseQualities(record, 25, 49, lowBaseQual);

        read = createRead(codingRegion, record, false, false);
        assertNotNull(read);

        assertEquals(25, read.trimmedBases());
        assertEquals(100, read.PositionStart);
        assertEquals(124, read.PositionEnd);
        assertEquals(0, read.ReadIndexStart);
        assertEquals(24, read.ReadIndexEnd);

        // repeat on the negative strand
        record.setReadNegativeStrandFlag(true);
        setBaseQualities(record, 0, 49, DEFAULT_MIN_BASE_QUAL); // reset

        setBaseQualities(record, 0, 9, lowBaseQual);

        codingRegion = new BaseRegion(100, 200);

        read = createRead(codingRegion, record, false, false);
        assertNotNull(read);

        assertEquals(10, read.trimmedBases());
        assertEquals(110, read.PositionStart);
        assertEquals(149, read.PositionEnd);
        assertEquals(10, read.ReadIndexStart);
        assertEquals(49, read.ReadIndexEnd);

        // also with soft-clipped bases being trimmed
        record = SamRecordTestUtils.createSamRecord(
                TEST_READ_ID, CHR_1, 120, TEST_READ_BASES.substring(0, 50), "20S30M", CHR_1, 300,
                false, false, null);

        record.setReadNegativeStrandFlag(true);
        setBaseQualities(record, 0, 9, lowBaseQual);

        read = createRead(codingRegion, record, true, false);
        assertNotNull(read);

        assertEquals(10, read.trimmedBases());
        assertEquals(110, read.PositionStart);
        assertEquals(149, read.PositionEnd);
        assertEquals(10, read.ReadIndexStart);
        assertEquals(49, read.ReadIndexEnd);
        assertEquals(10, read.SoftClippedStart);
        assertEquals(0, read.SoftClippedEnd);

        setBaseQualities(record, 0, 29, lowBaseQual);

        read = createRead(codingRegion, record, true, false);
        assertNotNull(read);

        assertEquals(30, read.trimmedBases());
        assertEquals(130, read.PositionStart);
        assertEquals(149, read.PositionEnd);
        assertEquals(30, read.ReadIndexStart);
        assertEquals(49, read.ReadIndexEnd);
        assertEquals(0, read.SoftClippedStart);
        assertEquals(0, read.SoftClippedEnd);
    }

    private static void setBaseQualities(final SAMRecord record, int rangeStart, int rangeEnd, byte baseQual)
    {
        for(int i = rangeStart; i <= rangeEnd; ++i)
        {
            record.getBaseQualities()[i] = baseQual;
        }
    }

    @Test
    public void testAdapterSoftClippedFragments()
    {
        BaseRegion codingRegion = new BaseRegion(100, 1000);

        // firstly a record with soft-clips at both ends which will be used
        SAMRecord record = SamRecordTestUtils.createSamRecord(
                TEST_READ_ID, CHR_1, 200, TEST_READ_BASES.substring(0, 80), "10S60M10S", CHR_1, 300,
                false, false, null);

        record.setInferredInsertSize(200);

        Read read = createRead(codingRegion, record, true, true);
        assertNotNull(read);

        assertEquals(190, read.PositionStart);
        assertEquals(269, read.PositionEnd);
        assertEquals(10, read.SoftClippedStart);
        assertEquals(10, read.SoftClippedEnd);

        // now restricted at the lower 3' end
        record.setInferredInsertSize(70);
        record.setReadNegativeStrandFlag(true);

        read = createRead(codingRegion, record, true, true);
        assertEquals(200, read.PositionStart);
        assertEquals(269, read.PositionEnd);
        assertEquals(0, read.SoftClippedStart);
        assertEquals(10, read.SoftClippedEnd);

        record.setReadNegativeStrandFlag(false);

        read = createRead(codingRegion, record, true, true);
        assertEquals(190, read.PositionStart);
        assertEquals(259, read.PositionEnd);
        assertEquals(10, read.SoftClippedStart);
        assertEquals(0, read.SoftClippedEnd);
    }

    @Test
    public void testFragmentAddNucleotide()
    {
        List<Integer> indices = Lists.newArrayList(1, 2, 3, 6, 7, 8);
        List<Byte> qualities = Lists.newArrayList((byte) 37, (byte) 25, (byte) 37, (byte) 37, (byte) 37, (byte) 25);
        List<String> nucleotides = Lists.newArrayList("A", "G", "T", "C", "A", "G");

        Fragment fragment = new Fragment(createReadRecord("01"), HLA_A, Sets.newHashSet(HLA_A), indices, qualities, nucleotides);

        fragment.removeLowQualBases();

        assertFalse(fragment.containsNucleotideLocus(2));

        fragment.addNucleotide(5, "G", (byte) 30);

        assertTrue(fragment.containsNucleotideLocus(5));

        fragment.addNucleotide(9, "T", (byte) 30);
        assertTrue(fragment.containsNucleotideLocus(9));
        assertTrue(fragment.validate());

        fragment.addNucleotide(0, "A", (byte) 37);
        assertTrue(fragment.containsNucleotideLocus(0));
        assertTrue(fragment.validate());
    }

    @Test
    public void testFragmentMerge()
    {
        String readId = "01";
        Read read = createReadRecord(readId);
        Fragment frag1 = new Fragment(
                read, HLA_A, Sets.newHashSet(HLA_A),
                Lists.newArrayList(1), Lists.newArrayList((byte) 30), Lists.newArrayList("A"));

        Fragment frag2 = new Fragment(
                read, HLA_B, Sets.newHashSet(HLA_B),
                Lists.newArrayList(1), Lists.newArrayList((byte) 30), Lists.newArrayList("A"));

        Fragment mergedFrag = mergeFragments(frag1, frag2);
        assertTrue(frag1.validate());
        assertEquals(2, mergedFrag.genes().size());
        assertEquals(1, mergedFrag.nucleotidesByLoci().size());
        assertEquals(1, mergedFrag.minNucleotideLocus());
        assertEquals(1, mergedFrag.nucleotidesByLoci().size());
        assertEquals(1, mergedFrag.nucleotidesByLoci().size());

        frag2 = new Fragment(
                read, HLA_A, Sets.newHashSet(HLA_A),
                Lists.newArrayList(0, 1, 2, 3),
                Lists.newArrayList((byte) 30, (byte) 30, (byte) 30, (byte) 30),
                Lists.newArrayList("A", "A", "A", "A"));

        mergedFrag = mergeFragments(frag1, frag2);
        assertTrue(frag1.validate());
        assertEquals(2, mergedFrag.genes().size());
        assertEquals(4, mergedFrag.nucleotidesByLoci().size());
        assertEquals(Integer.valueOf(0), Lists.newArrayList(mergedFrag.nucleotidesByLoci().keySet()).get(0));
        assertEquals(Integer.valueOf(1), Lists.newArrayList(mergedFrag.nucleotidesByLoci().keySet()).get(1));
        assertEquals(4, mergedFrag.nucleotidesByLoci().size());
        assertEquals(4, mergedFrag.nucleotidesByLoci().size());

        frag2 = new Fragment(
                read, HLA_C, Sets.newHashSet(HLA_C),
                Lists.newArrayList(3, 4, 5),
                Lists.newArrayList((byte) 30, (byte) 30, (byte) 30),
                Lists.newArrayList("A", "A", "A"));

        mergedFrag = mergeFragments(frag1, frag2);
        ArrayList<Integer> nucleotideLoci = Lists.newArrayList(mergedFrag.nucleotidesByLoci().keySet());
        assertTrue(frag1.validate());
        assertEquals(3, mergedFrag.genes().size());
        assertEquals(6, mergedFrag.nucleotidesByLoci().size());
        assertEquals(Integer.valueOf(0), nucleotideLoci.get(0));
        assertEquals(Integer.valueOf(3), nucleotideLoci.get(3));
        assertEquals(Integer.valueOf(4), nucleotideLoci.get(4));
        assertEquals(Integer.valueOf(5), nucleotideLoci.get(5));
        assertEquals(6, mergedFrag.nucleotidesByLoci().size());
        assertEquals(6, mergedFrag.nucleotidesByLoci().size());
    }

    private static void assertRange(int expectedStart, int expectedEnd, List<Integer> victim)
    {
        assertEquals(expectedStart, victim.get(0).intValue());
        assertEquals(expectedEnd, victim.get(victim.size() - 1).intValue());
    }

    @Ignore
    @Test
    public void testLociLookup()
    {
        // test locus look-up for a fragment with typical number of nucleotides (~150 being a read within an exon, high value = 250)
        // relative search times: manual = 5, binary = 2, indexOf = 3
        int iterations = 1_000_000;

        int lociCount = 150;

        List<Integer> lociValues = Lists.newArrayListWithCapacity(lociCount);

        for(int i = 0; i < lociCount; ++i)
        {
            lociValues.add(i);
        }

        int[] searchLoci = {0, 10, 25, 50, 100, lociCount - 1};
        int searchValueIndex = 0;

        long startTime = System.nanoTime();

        long sampleTime = System.nanoTime() - startTime;
        double sampleTimeMillis = sampleTime / NANOS_IN_MILLISECOND;

        LL_LOGGER.info(format("manual search time: %.6fms", sampleTimeMillis));

        startTime = System.nanoTime();

        sampleTime = System.nanoTime() - startTime;
        sampleTimeMillis = sampleTime / NANOS_IN_MILLISECOND;

        LL_LOGGER.info(format("binary search time: %.6fms", sampleTimeMillis));

        startTime = System.nanoTime();

        for(int i = 0; i < iterations; ++i)
        {
            lociValues.indexOf(searchLoci[searchValueIndex++]);

            if(searchValueIndex >= searchLoci.length)
                searchValueIndex = 0;
        }

        sampleTime = System.nanoTime() - startTime;
        sampleTimeMillis = sampleTime / NANOS_IN_MILLISECOND;

        LL_LOGGER.info(format("indexOf search time: %.6fms", sampleTimeMillis));
    }
}
