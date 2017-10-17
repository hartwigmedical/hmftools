package com.hartwig.hmftools.common.purple.copynumber;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleCopyNumberTest {

    @Test
    public void testNegativeCopyNumber() {
        PurpleCopyNumber copyNumber = PurpleDatamodelTest.createCopyNumber("1", 1, 100, -12).build();
        assertEquals(0, copyNumber.value());
    }

    @Test
    public void testNegativeBaf() {
        assertEquals("AAAAAAA", createBAFString(-0.094, 7));
    }

    @Test
    public void testDescriptiveBAFBoundary() {
        assertEquals("AA", createBAFString(1.3, 2));
        assertEquals("", createBAFString(0, -12));
    }

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

        assertEquals("A[12]x", createBAFString(0, 12));
    }

    @NotNull
    private String createBAFString(double actualBaf, double copyNumber) {
        return PurpleDatamodelTest.createCopyNumber("1", 1, 100, copyNumber)
                .averageObservedBAF(actualBaf)
                .averageActualBAF(actualBaf)
                .build()
                .descriptiveBAF();
    }
}
