package com.hartwig.hmftools.common.bam;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import org.junit.Test;

public class SupplementaryReadDataTest
{
    private static final SupplementaryReadData TEST_SUPP_DATA1 = new SupplementaryReadData(
            "chr20", 31055958, SupplementaryReadData.SUPP_POS_STRAND, "20S14M5I81M31S", 1, 8);

    private static final SupplementaryReadData TEST_SUPP_DATA2 = new SupplementaryReadData(
            "chrUn_GL000216v2", 13986, SupplementaryReadData.SUPP_NEG_STRAND, "116S35M", 0, 0);

    @Test
    public void testFromAlignmentEmptyAlignmentReturnsNull()
    {
        assertNull(SupplementaryReadData.fromAlignment(""));
    }

    @Test
    public void testFromAlignmentActualAlignment()
    {
        final SupplementaryReadData actual = SupplementaryReadData.fromAlignment(TEST_SUPP_DATA1.asSamTag());

        assertNotNull(actual);
        assertEquals(TEST_SUPP_DATA1, actual);
    }

    @Test
    public void testFromAlignmentAlternativeDelimiter()
    {
        final SupplementaryReadData actual = SupplementaryReadData.fromAlignment(TEST_SUPP_DATA1.asDelimStr(), ";");

        assertNotNull(actual);
        assertEquals(TEST_SUPP_DATA1, actual);
    }

    @Test
    public void testFromEmtpyReturnsNull()
    {
        assertNull(SupplementaryReadData.extractAlignments(""));
    }

    @Test
    public void testFromSingleAlignmentWithEndingDelim()
    {
        final List<SupplementaryReadData> actual = SupplementaryReadData.extractAlignments(TEST_SUPP_DATA1.asSamTag() + ';');

        assertNotNull(actual);
        assertEquals(1, actual.size());
        assertEquals(TEST_SUPP_DATA1, actual.get(0));
    }

    @Test
    public void testFromSingleAlignmentWithoutEndingDelim()
    {
        final List<SupplementaryReadData> actual = SupplementaryReadData.extractAlignments(TEST_SUPP_DATA1.asSamTag());

        assertNotNull(actual);
        assertEquals(1, actual.size());
        assertEquals(TEST_SUPP_DATA1, actual.get(0));
    }

    @Test
    public void testFromSingleAlignmentAlternativeDelimiter()
    {
        final List<SupplementaryReadData> actual = SupplementaryReadData.extractAlignments(TEST_SUPP_DATA1.asDelimStr());

        assertNotNull(actual);
        assertEquals(1, actual.size());
        assertEquals(TEST_SUPP_DATA1, actual.get(0));
    }

    @Test
    public void testFromTwoAlignmentsWithEndingDelim()
    {
        final List<SupplementaryReadData> actual =
                SupplementaryReadData.extractAlignments(TEST_SUPP_DATA1.asSamTag() + ';' + TEST_SUPP_DATA2.asSamTag() + ';');

        assertNotNull(actual);
        assertEquals(2, actual.size());
        assertEquals(TEST_SUPP_DATA1, actual.get(0));
        assertEquals(TEST_SUPP_DATA2, actual.get(1));
    }

    @Test
    public void testFromTwoAlignmentsWithoutEndingDelim()
    {
        final List<SupplementaryReadData> actual =
                SupplementaryReadData.extractAlignments(TEST_SUPP_DATA1.asSamTag() + ';' + TEST_SUPP_DATA2.asSamTag());

        assertNotNull(actual);
        assertEquals(2, actual.size());
        assertEquals(TEST_SUPP_DATA1, actual.get(0));
        assertEquals(TEST_SUPP_DATA2, actual.get(1));
    }

    @Test
    public void testFromTwoAlignmentEmptyFirstReturnsNull()
    {
        assertNull(SupplementaryReadData.extractAlignments(";" + TEST_SUPP_DATA1.asSamTag()));
    }

    @Test
    public void testFromTwoAlignmentEmptyLastReturnsNull()
    {
        assertNull(SupplementaryReadData.extractAlignments(TEST_SUPP_DATA1.asSamTag() + ";;"));
    }

    @Test
    public void testFirstAlignmentFromEmtpyReturnsNull()
    {
        assertNull(SupplementaryReadData.extractAlignment(""));
    }

    @Test
    public void testFirstAlignmentFromSingleAlignmentWithEndingDelim()
    {
        final SupplementaryReadData actual = SupplementaryReadData.extractAlignment(TEST_SUPP_DATA1.asSamTag() + ';');

        assertNotNull(actual);
        assertEquals(TEST_SUPP_DATA1, actual);
    }

    @Test
    public void testFirstAlignmentFromSingleAlignmentWithoutEndingDelim()
    {
        final SupplementaryReadData actual = SupplementaryReadData.extractAlignment(TEST_SUPP_DATA1.asSamTag());

        assertNotNull(actual);
        assertEquals(TEST_SUPP_DATA1, actual);
    }

    @Test
    public void testFirstAlignmentFromSingleAlignmentAlternativeDelimiter()
    {
        final SupplementaryReadData actual = SupplementaryReadData.extractAlignment(TEST_SUPP_DATA1.asDelimStr());

        assertNotNull(actual);
        assertEquals(TEST_SUPP_DATA1, actual);
    }

    @Test
    public void testFirstAlignmentFromTwoAlignmentsWithEndingDelim()
    {
        final SupplementaryReadData actual =
                SupplementaryReadData.extractAlignment(TEST_SUPP_DATA1.asSamTag() + ';' + TEST_SUPP_DATA2.asSamTag() + ';');

        assertNotNull(actual);
        assertEquals(TEST_SUPP_DATA1, actual);
    }

    @Test
    public void testFirstAlignmentFromTwoAlignmentsWithoutEndingDelim()
    {
        final SupplementaryReadData actual =
                SupplementaryReadData.extractAlignment(TEST_SUPP_DATA1.asSamTag() + ';' + TEST_SUPP_DATA2.asSamTag());

        assertNotNull(actual);
        assertEquals(TEST_SUPP_DATA1, actual);
    }

    @Test
    public void testFirstAlignmentFromTwoAlignmentEmptyFirst_ReturnsNull()
    {
        assertNull(SupplementaryReadData.extractAlignment(";" + TEST_SUPP_DATA1.asSamTag()));
    }

    @Test
    public void testFirstAlignmentFromTwoAlignmentEmptyLast()
    {
        final SupplementaryReadData actual = SupplementaryReadData.extractAlignment(TEST_SUPP_DATA1.asSamTag() + ";;");

        assertNotNull(actual);
        assertEquals(TEST_SUPP_DATA1, actual);
    }

    @Test
    public void testAlignmentCountTwoAlignmentsEndingWithDelim()
    {
        final String suppData = TEST_SUPP_DATA1.asSamTag() + ';' + TEST_SUPP_DATA2.asSamTag() + ';';

        assertEquals(2, SupplementaryReadData.alignmentCount(suppData));
    }

    @Test
    public void testAlignmentCountTwoAlignmentsEndingWithoutDelim()
    {
        final String suppData = TEST_SUPP_DATA1.asSamTag() + ';' + TEST_SUPP_DATA2.asSamTag();

        assertEquals(2, SupplementaryReadData.alignmentCount(suppData));
    }
}
