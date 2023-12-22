package com.hartwig.hmftools.svassembly.assembly;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;

import com.hartwig.hmftools.svassembly.TestUtils;
import com.hartwig.hmftools.svassembly.models.Record;

import org.junit.Test;

public class NodeFolderTest
{
    @Test
    public void canFoldPaths()
    {
        final List<Record> sampleData = List.of(
                TestUtils.createSAMRecord("AATAA"),
                TestUtils.createSAMRecord("AAAAA"),
                TestUtils.createSAMRecord("AAAAA"),
                TestUtils.createSAMRecord("AAAAA"),
                TestUtils.createSAMRecord("AAAAA"),
                TestUtils.createSAMRecord("AAAAA")
        );

        final HeadNode merged = sampleData.stream().map(record -> HeadNode.create(TestUtils.config(), record, 1, true))
                .reduce(HeadNode::combine)
                .orElseThrow();

        new NodeFolder(TestUtils.config()).foldPaths(merged);

        assertThat(merged.flatten()).containsOnly("AAAAA", "AATAA");
        assertThat(merged.Base).isEqualTo('S');
        //noinspection DataFlowIssue
        assertThat(merged.nextA.nextA.nextA.nextA).isSameAs(merged.nextA.nextA.nextT.nextA);
    }

    @Test
    public void canFoldPathsBranching()
    {
        final List<Record> sampleData = List.of(
                TestUtils.createSAMRecord("AATAAAGC"),
                TestUtils.createSAMRecord("AAAAAAGG"),
                TestUtils.createSAMRecord("AAAAAAGG"),
                TestUtils.createSAMRecord("AAAAAAGG"),
                TestUtils.createSAMRecord("AAAAAAGG"),
                TestUtils.createSAMRecord("AAAAAAGG")
        );
        sampleData.get(0).getBaseQuality()["AATAAAGC".length() - 1] = 11;

        final HeadNode merged = sampleData.stream().map(record -> HeadNode.create(TestUtils.config(), record, 1, true))
                .reduce(HeadNode::combine)
                .orElseThrow();

        new NodeFolder(TestUtils.config()).foldPaths(merged);

        assertThat(merged.Base).isEqualTo('S');
        assertThat(merged.flatten()).containsExactlyInAnyOrder("AAAAAAGG", "AATAAAGG", "AAAAAAGC", "AATAAAGC");
        //noinspection DataFlowIssue
        assertThat(merged.nextA.nextA.nextA.nextA).isSameAs(merged.nextA.nextA.nextT.nextA);
    }
}
