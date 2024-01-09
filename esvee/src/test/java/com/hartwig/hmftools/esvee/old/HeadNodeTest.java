package com.hartwig.hmftools.esvee.old;

import static com.hartwig.hmftools.esvee.TestUtils.createSAMRecord;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Objects;

import com.hartwig.hmftools.esvee.TestUtils;
import com.hartwig.hmftools.esvee.old.HeadNode;
import com.hartwig.hmftools.esvee.read.Read;

import org.junit.Test;

public class HeadNodeTest
{
    private static HeadNode create(final String sequence, final int position, final boolean isForwards)
    {
        return HeadNode.create(TestUtils.createSAMRecord(sequence), position, isForwards);
    }

    private static HeadNode create(final Read read, final int position, final boolean isForwards)
    {
        return Objects.requireNonNull(HeadNode.create(read, position, isForwards));
    }

    @Test
    public void createGraphsFromSequences()
    {
        assertEquals("AATAA", create("AATAA", 1, true).flatten().get(0));

        assertEquals("AG", create("AAAAG", 4, true).flatten().get(0));

        assertEquals("A", create("AATAA", 1, false).flatten().get(0));

        assertEquals("AAAAT", create("TAAAA", 5, false).flatten().get(0));
    }

    @Test
    public void testBasicNodeRoutines()
    {
        final int junctionPos = 10;

        final List<Read> records = List.of(
                createSAMRecord("AACCGG", 10),
                createSAMRecord("AATCGG", 10),
                createSAMRecord("ACCCGG", 10));

        HeadNode merged = records.stream().map(record -> HeadNode.create(record, junctionPos, true))
                .reduce(HeadNode::combine)
                .orElseThrow();

        // assertEquals(3, merged.junctionCount());

        // merged.trimShortPaths(10);

        List<String> flattened = merged.flatten();

        assertEquals(3, flattened.size());
        assertTrue(flattened.contains("AACCGG"));
        assertTrue(flattened.contains("AATCGG"));
        assertTrue(flattened.contains("ACCCGG"));

    }

        @Test
    public void clipsShortPaths()
    {
        final int junctionPos = 10;
        final List<Read> records = List.of(
                createSAMRecord("AAAAAAAAAAAA", 10),
                createSAMRecord("TTTTTTTTTTTT", 10),
                createSAMRecord("CCCCCCCCC", 10),
                createSAMRecord("GGGGGGGGGG", 10),
                createSAMRecord("ATATAT", 10),
                createSAMRecord("ATATACGGGGGG", junctionPos), // branching of the one above
                createSAMRecord("CGCGCGCGCGCG", 10)
        );

        HeadNode merged = records.stream().map(record -> HeadNode.create(record, junctionPos, true))
                .reduce(HeadNode::combine)
                .orElseThrow();

        merged.trimShortPaths(10);

        List<String> flattened = merged.flatten();

        assertFalse(flattened.contains("ATATAT"));

        assertTrue(flattened.contains("AAAAAAAAAAAA"));
        assertTrue(flattened.contains("TTTTTTTTTTTT"));
        assertTrue(flattened.contains("GGGGGGGGGG"));
        assertTrue(flattened.contains("ATATACGGGGGG"));
        assertTrue(flattened.contains("CGCGCGCGCGCG"));
    }

    @Test
    public void ignoresConfidentIncorrectBasesAtStart()
    {
        int junctionPos = 10;

        String mainSequence = "ATCGATCGATCG";
        String altSequence = "TTCGATCGATCG";

        List<Read> records = List.of(
                TestUtils.createSAMRecord(mainSequence, junctionPos),
                TestUtils.createSAMRecord(altSequence, junctionPos),
                TestUtils.createSAMRecord(mainSequence, junctionPos),
                TestUtils.createSAMRecord(mainSequence, junctionPos),
                TestUtils.createSAMRecord(mainSequence, junctionPos),
                TestUtils.createSAMRecord(mainSequence, junctionPos));

        HeadNode merged = records.stream().map(record -> HeadNode.create(record, junctionPos, true))
                .reduce(HeadNode::combine)
                .orElseThrow();

        List<String> flattened = merged.flatten();

        assertEquals(2, flattened.size());
        assertTrue(flattened.contains(mainSequence));
        assertTrue(flattened.contains(altSequence));
    }

    /* CHASHA FIXME


    @Test
    public void createsAlignmentsContainingInsertsCorrectlyForwards() {
        //                       Junction v
        //                       SSSMMMMIMMMM
        // Reference / Alignment 7890123.4567
        final String sequence = "AAATTTCCCGGG";
        final Record record = TestUtils.createSAMRecord(sequence, 10);
        record.setCigar("3S4M1I4M");

        assertTrue(create(record, 15, true).flatten()).containsExactly("GGG");
    }

    @Test
    public void createsAlignmentsContainingInsertsCorrectlyReverse() {
        //                       Junction v
        //                       SSSMMMMIMMMM
        // Reference / Alignment 7890123.4567
        final String sequence = "AAATTTCCCGGG";
        final Record record = TestUtils.createSAMRecord(sequence, 10);
        record.setCigar("3S4M1I4M");

        assertTrue(create(record, 15, false).flatten()).containsExactly("GCCCTTTAAA");
    }

    // the following are for AssemblyExtension

    @Test
    public void canAttachCompletelyOverlapping()
    {
        final HeadNode existing = create("AAATTTCCCGGGAATTCCGGATCG", 1, true);
        final HeadNode incoming = create("AAATTTCCCGGGAATTCCGGATCGAA", 1, true);

        assertTrue(existing.attach(incoming, 10, 1, 0, Integer.MAX_VALUE)).isTrue();
        assertTrue(existing.flatten()).hasSize(1);
        assertTrue(existing.flatten()).containsExactly("AAATTTCCCGGGAATTCCGGATCGAA");
    }

    @Test
    public void canAttachSingleCandidatePoint()
    {
        final HeadNode existing = create("AAATTTCCCGGGAATTCCGGATCG", 1, true);
        final HeadNode incoming = create("AATTTCCCGGGAATTCCGGATCGAA", 1, true);

        assertTrue(existing.attach(incoming, 10, 1, 0, Integer.MAX_VALUE)).isTrue();
        assertTrue(existing.flatten()).hasSize(1);
        assertTrue(existing.flatten()).containsExactly("AAATTTCCCGGGAATTCCGGATCGAA");
    }

    @Test
    public void canAttachSingleCandidatePointMinorError()
    {
        final HeadNode existing = create("AAATTTCCCGGGAATTCCGGATCG", 1, true);
        final HeadNode incoming = create("AATTTCCCGGGATTTCCGGATCGAA", 1, true);

        assertTrue(existing.attach(incoming, 10, 1, 0, Integer.MAX_VALUE)).isTrue();
        assertTrue(existing.flatten()).hasSize(2);
        assertTrue(existing.flatten()).containsExactlyInAnyOrder("AAATTTCCCGGGAATTCCGGATCGAA", "AAATTTCCCGGGATTTCCGGATCGAA");
    }

    @Test
    public void canAttachSingleCandidatePointLowQualErrors()
    {
        final HeadNode existing = create("AAATTTCCCGGGAATTCCGGATCG", 1, true);

        final Record incomingRecord = createSAMRecord("AAAAACCCGGGAATTCCGGATCGAA");
        final byte[] baseQualities = incomingRecord.getBaseQuality();
        for(int i = 2; i < 5; i++)
            baseQualities[i] = 11;

        final HeadNode incoming = create(incomingRecord, 1, true);

        assertTrue(existing.attach(incoming, 10, 1, 0, Integer.MAX_VALUE)).isTrue();

        new NodeFolder(TestUtils.config()).foldPaths(existing);
        existing.pruneNodes();

        assertTrue(existing.flatten()).hasSize(1);
        assertTrue(existing.flatten()).containsExactly("AAATTTCCCGGGAATTCCGGATCGAA");
    }

    @Ignore
    @Test
    public void canAttachMultipleCandidatePoints()
    {
        final HeadNode existing = create("AAAATTTCCCAAAATTTCCC", 1, true);
        final HeadNode incoming = create("AAAATTTCCCAAAAC", 1, true);

        assertTrue(existing.attach(incoming, 10, 1, 0, Integer.MAX_VALUE)).isTrue();
        assertTrue(existing.flatten()).hasSize(2);
        assertTrue(existing.flatten()).containsExactlyInAnyOrder("AAAATTTCCCAAAAC", "AAAATTTCCCAAAATTTCCCAAAAC");
    }

    @Test
    public void doesNotAttachAmbiguous()
    {
        final HeadNode existing = create("ATATATATATATATATATAT", 1, true);
        final HeadNode incoming = create("ATATATATATAT", 1, true);

        assertTrue(existing.attach(incoming, 10, 1, 0, Integer.MAX_VALUE)).isFalse();
    }

    @Test
    public void doesNotAttachUnrelated()
    {
        final HeadNode existing = create("AAATTTCCCGGGAAATTTCCCGGG", 1, true);
        final HeadNode incoming = create("AAACCCTTTGGGAAACCCTTT", 1, true);

        assertTrue(existing.attach(incoming, 10, 1, 0, Integer.MAX_VALUE)).isFalse();
    }

    @Test
    public void doesNotAttachIfMoreThanOneError()
    {
        final HeadNode existing = create("AAATTTCCCGGGAAATTTCCCGGG", 1, true);
        final HeadNode incoming = create("AAATATCTCGGGAAATTTCCCGGGAA", 1, true);

        assertTrue(existing.attach(incoming, 10, 1, 0, Integer.MAX_VALUE)).isFalse();
    }

    @Test
    public void doesNotAttachIfLessThanMinOverlap()
    {
        final HeadNode existing = create("AAATTTCCCGGGAAATTTCCCGGG", 1, true);
        final HeadNode incoming = create("TTTCCCGGGTT", 1, true);

        assertTrue(existing.attach(incoming, 10, 1, 0, Integer.MAX_VALUE)).isFalse();
    }
    */
}
