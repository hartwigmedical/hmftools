package com.hartwig.hmftools.common.purple.copynumber;

import static org.junit.Assert.assertEquals;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleCopyNumberTest {

    @Test
    public void testDescriptiveBAF() {
        assertEquals("AA", createBAFString(0, 2));
        assertEquals("AB", createBAFString(0.5, 2));
        assertEquals("AB", createBAFString(0.7, 2));
        assertEquals("AA", createBAFString(0.8, 2));
        assertEquals("AA", createBAFString(1, 2));

        assertEquals("AAA", createBAFString(0, 3));
        assertEquals("AAB", createBAFString(0.5, 3));
        assertEquals("AAB", createBAFString(0.83, 3));
        assertEquals("AAA", createBAFString(0.84, 3));
        assertEquals("AAA", createBAFString(1, 3));

        assertEquals("AAAA", createBAFString(0, 4));
        assertEquals("AABB", createBAFString(0.5, 4));
        assertEquals("AABB", createBAFString(0.624, 4));
        assertEquals("AAAB", createBAFString(0.625, 4));
        assertEquals("AAAB", createBAFString(0.874, 4));
        assertEquals("AAAA", createBAFString(0.875, 4));
        assertEquals("AAAA", createBAFString(0.9, 4));
    }

    @NotNull
    private String createBAFString(double actualBaf, double copyNumber) {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome("1")
                .start(1)
                .end(100)
                .averageTumorCopyNumber(copyNumber)
                .bafCount(0)
                .averageObservedBAF(actualBaf)
                .averageActualBAF(actualBaf)
                .build().descriptiveBAF();
    }
}
