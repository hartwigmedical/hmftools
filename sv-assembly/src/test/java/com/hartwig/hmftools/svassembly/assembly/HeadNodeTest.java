package com.hartwig.hmftools.svassembly.assembly;

import static com.hartwig.hmftools.svassembly.TestUtils.createSAMRecord;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;
import java.util.Objects;

import com.hartwig.hmftools.svassembly.TestUtils;
import com.hartwig.hmftools.svassembly.models.Record;

import org.junit.Ignore;
import org.junit.Test;

public class HeadNodeTest
{
    private static HeadNode create(final String sequence, final int position, final boolean isForwards)
    {
        return HeadNode.create(TestUtils.config(), createSAMRecord(sequence), position, isForwards);
    }

    private static HeadNode create(final Record record, final int position, final boolean isForwards)
    {
        return Objects.requireNonNull(HeadNode.create(TestUtils.config(), record, position, isForwards));
    }

    @Test
    public void createGraphsFromSequencesForwards()
    {
        assertThat(create("AATAA", 1, true).flatten().get(0)).isEqualTo("AATAA");
        assertThat(create("TAAAA", 1, true).flatten().get(0))
                .isEqualTo("TAAAA");
        assertThat(create("AAAAG", 1, true).flatten().get(0))
                .isEqualTo("AAAAG");

        assertThat(create("AAAAG", 4, true).flatten().get(0))
                .isEqualTo("AG");
    }

    @Test
    public void createGraphsFromSequencesBackwards()
    {
        assertThat(create("AATAA", 1, false).flatten().get(0))
                .isEqualTo("A");
        assertThat(create("TAAAA", 5, false).flatten().get(0))
                .isEqualTo("AAAAT");
        assertThat(create("AAAAG", 5, false).flatten().get(0))
                .isEqualTo("GAAAA");

        assertThat(create("AAAAG", 4, false).flatten().get(0))
                .isEqualTo("AAAA");
        assertThat(create("TAAAG", 4, false).flatten().get(0))
                .isEqualTo("AAAT");
    }

    @Test
    public void clipsShortPaths()
    {
        final int junctionPos = 10;
        final List<Record> records = List.of(
                createSAMRecord("AAAAAAAAAAAA", 10),
                createSAMRecord("TTTTTTTTTTTT", 10),
                createSAMRecord("CCCCCCCCC", 10),
                createSAMRecord("GGGGGGGGGG", 10),
                createSAMRecord("ATATAT", 10),
                createSAMRecord("CGCGCGCGCGCG", 10)
        );

        final HeadNode merged = records.stream().map(record -> HeadNode.create(TestUtils.config(), record, junctionPos, true))
                .reduce(HeadNode::combine)
                .orElseThrow();

        merged.trimShortPaths(10);
        final List<String> flattened = merged.flatten();
        assertThat(flattened).doesNotContain("ATATAT");
        assertThat(flattened).containsExactlyInAnyOrder(
                "AAAAAAAAAAAA",
                "TTTTTTTTTTTT",
                "GGGGGGGGGG",
                "CGCGCGCGCGCG"
        );
    }

    @Test
    public void clipsShortPathsBranching()
    {
        final int junctionPos = 10;
        final List<Record> records = List.of(
                createSAMRecord("AAAAAAAAAAAA", junctionPos),
                createSAMRecord("TTTTTTTTTTTT", junctionPos),
                createSAMRecord("CCCCCCCCCCCC", junctionPos),
                createSAMRecord("ATATACGGGGGG", junctionPos),
                createSAMRecord("ATATAT", junctionPos),
                createSAMRecord("CGCGCGCGCGCG", junctionPos));

        final HeadNode merged = records.stream().map(record -> HeadNode.create(TestUtils.config(), record, junctionPos, true))
                .reduce(HeadNode::combine)
                .orElseThrow();

        merged.trimShortPaths(10);
        final List<String> flattened = merged.flatten();
        assertThat(flattened).doesNotContain("ATATAT");
        assertThat(flattened).containsExactlyInAnyOrder(
                "AAAAAAAAAAAA",
                "TTTTTTTTTTTT",
                "CCCCCCCCCCCC",
                "ATATACGGGGGG",
                "CGCGCGCGCGCG"
        );
    }

    @Test
    public void ignoresConfidentIncorrectBasesAtStart()
    {
        final int junctionPos = 10;
        final List<Record> records = List.of(
                TestUtils.createSAMRecord("ATCGATCGATCG", junctionPos),
                TestUtils.createSAMRecord("TTCGATCGATCG", junctionPos),
                TestUtils.createSAMRecord("ATCGATCGATCG", junctionPos),
                TestUtils.createSAMRecord("ATCGATCGATCG", junctionPos),
                TestUtils.createSAMRecord("ATCGATCGATCG", junctionPos),
                TestUtils.createSAMRecord("ATCGATCGATCG", junctionPos));

        final HeadNode merged = records.stream().map(record -> HeadNode.create(TestUtils.config(), record, junctionPos, true))
                .reduce(HeadNode::combine)
                .orElseThrow();
        // TODO: Process this
    }

    @Test
    public void createsAlignmentsContainingInsertsCorrectlyForwards() {
        //                       Junction v
        //                       SSSMMMMIMMMM
        // Reference / Alignment 7890123.4567
        final String sequence = "AAATTTCCCGGG";
        final Record record = TestUtils.createSAMRecord(sequence, 10);
        record.setCigar("3S4M1I4M");

        assertThat(create(record, 15, true).flatten()).containsExactly("GGG");
    }

    @Test
    public void createsAlignmentsContainingInsertsCorrectlyReverse() {
        //                       Junction v
        //                       SSSMMMMIMMMM
        // Reference / Alignment 7890123.4567
        final String sequence = "AAATTTCCCGGG";
        final Record record = TestUtils.createSAMRecord(sequence, 10);
        record.setCigar("3S4M1I4M");

        assertThat(create(record, 15, false).flatten()).containsExactly("GCCCTTTAAA");
    }

    @Test
    public void canAttachCompletelyOverlapping()
    {
        final HeadNode existing = create("AAATTTCCCGGGAATTCCGGATCG", 1, true);
        final HeadNode incoming = create("AAATTTCCCGGGAATTCCGGATCGAA", 1, true);

        assertThat(existing.attach(incoming, 10, 1, 0, Integer.MAX_VALUE)).isTrue();
        assertThat(existing.flatten()).hasSize(1);
        assertThat(existing.flatten()).containsExactly("AAATTTCCCGGGAATTCCGGATCGAA");
    }

    @Test
    public void canAttachSingleCandidatePoint()
    {
        final HeadNode existing = create("AAATTTCCCGGGAATTCCGGATCG", 1, true);
        final HeadNode incoming = create("AATTTCCCGGGAATTCCGGATCGAA", 1, true);

        assertThat(existing.attach(incoming, 10, 1, 0, Integer.MAX_VALUE)).isTrue();
        assertThat(existing.flatten()).hasSize(1);
        assertThat(existing.flatten()).containsExactly("AAATTTCCCGGGAATTCCGGATCGAA");
    }

    @Test
    public void canAttachSingleCandidatePointMinorError()
    {
        final HeadNode existing = create("AAATTTCCCGGGAATTCCGGATCG", 1, true);
        final HeadNode incoming = create("AATTTCCCGGGATTTCCGGATCGAA", 1, true);

        assertThat(existing.attach(incoming, 10, 1, 0, Integer.MAX_VALUE)).isTrue();
        assertThat(existing.flatten()).hasSize(2);
        assertThat(existing.flatten()).containsExactlyInAnyOrder("AAATTTCCCGGGAATTCCGGATCGAA", "AAATTTCCCGGGATTTCCGGATCGAA");
    }

    @Test
    public void canAttachSingleCandidatePointLowQualErrors()
    {
        final HeadNode existing = create("AAATTTCCCGGGAATTCCGGATCG", 1, true);

        final Record incomingRecord = createSAMRecord("AAAAACCCGGGAATTCCGGATCGAA");
        final byte[] baseQualities = incomingRecord.getBaseQuality();
        for (int i = 2; i < 5; i++)
            baseQualities[i] = 11;

        final HeadNode incoming = create(incomingRecord, 1, true);

        assertThat(existing.attach(incoming, 10, 1, 0, Integer.MAX_VALUE)).isTrue();

        new NodeFolder(TestUtils.config()).foldPaths(existing);
        existing.pruneNodes();

        assertThat(existing.flatten()).hasSize(1);
        assertThat(existing.flatten()).containsExactly("AAATTTCCCGGGAATTCCGGATCGAA");
    }

    @Ignore
    @Test
    public void canAttachMultipleCandidatePoints()
    {
        final HeadNode existing = create("AAAATTTCCCAAAATTTCCC", 1, true);
        final HeadNode incoming = create("AAAATTTCCCAAAAC", 1, true);

        assertThat(existing.attach(incoming, 10, 1, 0, Integer.MAX_VALUE)).isTrue();
        assertThat(existing.flatten()).hasSize(2);
        assertThat(existing.flatten()).containsExactlyInAnyOrder("AAAATTTCCCAAAAC", "AAAATTTCCCAAAATTTCCCAAAAC");
    }

    @Test
    public void doesNotAttachAmbiguous()
    {
        final HeadNode existing = create("ATATATATATATATATATAT", 1, true);
        final HeadNode incoming = create("ATATATATATAT", 1, true);

        assertThat(existing.attach(incoming, 10, 1, 0, Integer.MAX_VALUE)).isFalse();
    }

    @Test
    public void doesNotAttachUnrelated()
    {
        final HeadNode existing = create("AAATTTCCCGGGAAATTTCCCGGG", 1, true);
        final HeadNode incoming = create("AAACCCTTTGGGAAACCCTTT", 1, true);

        assertThat(existing.attach(incoming, 10, 1, 0, Integer.MAX_VALUE)).isFalse();
    }

    @Test
    public void doesNotAttachIfMoreThanOneError()
    {
        final HeadNode existing = create("AAATTTCCCGGGAAATTTCCCGGG", 1, true);
        final HeadNode incoming = create("AAATATCTCGGGAAATTTCCCGGGAA", 1, true);

        assertThat(existing.attach(incoming, 10, 1, 0, Integer.MAX_VALUE)).isFalse();
    }

    @Test
    public void doesNotAttachIfLessThanMinOverlap()
    {
        final HeadNode existing = create("AAATTTCCCGGGAAATTTCCCGGG", 1, true);
        final HeadNode incoming = create("TTTCCCGGGTT", 1, true);

        assertThat(existing.attach(incoming, 10, 1, 0, Integer.MAX_VALUE)).isFalse();
    }
}
