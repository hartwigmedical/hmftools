package com.hartwig.hmftools.esvee.assembly;

public class NodeFolderTest
{
    /* CHASHA FIXME
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

        assertTrue(merged.flatten()).containsOnly("AAAAA", "AATAA");
        assertTrue(merged.Base).isEqualTo('S');
        //noinspection DataFlowIssue
        assertTrue(merged.nextA.nextA.nextA.nextA).isSameAs(merged.nextA.nextA.nextT.nextA);
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

        assertTrue(merged.Base).isEqualTo('S');
        assertTrue(merged.flatten()).containsExactlyInAnyOrder("AAAAAAGG", "AATAAAGG", "AAAAAAGC", "AATAAAGC");
        //noinspection DataFlowIssue
        assertTrue(merged.nextA.nextA.nextA.nextA).isSameAs(merged.nextA.nextA.nextT.nextA);
    }

     */
}
