package com.hartwig.hmftools.common.bam;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import org.junit.Test;

public class SupplementaryReadDataTest
{
    private static final SupplementaryReadData TEST_SUPP_DATA1 =
            new SupplementaryReadData("chr20", 31055958, SupplementaryReadData.SUPP_POS_STRAND, "20S14M5I81M31S", 1, 8);
    private static final SupplementaryReadData TEST_SUPP_DATA2 =
            new SupplementaryReadData("chrUn_GL000216v2", 13986, SupplementaryReadData.SUPP_NEG_STRAND, "116S35M", 0, 0);

    @Test
    public void testFromAlignment_EmptyAlignment_ReturnsNull()
    {
        assertNull(SupplementaryReadData.fromAlignment(""));
    }

    @Test
    public void testFromAlignment_ActualAlignment()
    {
        final SupplementaryReadData actual = SupplementaryReadData.fromAlignment(TEST_SUPP_DATA1.asSamTag());

        assertNotNull(actual);
        assertEquals(TEST_SUPP_DATA1, actual);
    }

    @Test
    public void testFromAlignment_AlternativeDelimiter()
    {
        final SupplementaryReadData actual = SupplementaryReadData.fromAlignment(TEST_SUPP_DATA1.asCsv(), ";");

        assertNotNull(actual);
        assertEquals(TEST_SUPP_DATA1, actual);
    }

    @Test
    public void testFrom_Emtpy_ReturnsNull()
    {
        assertNull(SupplementaryReadData.extractAlignments(""));
    }

    @Test
    public void testFrom_SingleAlignment_WithEndingDelim()
    {
        final List<SupplementaryReadData> actual = SupplementaryReadData.extractAlignments(TEST_SUPP_DATA1.asSamTag() + ';');

        assertNotNull(actual);
        assertEquals(1, actual.size());
        assertEquals(TEST_SUPP_DATA1, actual.get(0));
    }

    @Test
    public void testFrom_SingleAlignment_WithoutEndingDelim()
    {
        final List<SupplementaryReadData> actual = SupplementaryReadData.extractAlignments(TEST_SUPP_DATA1.asSamTag());

        assertNotNull(actual);
        assertEquals(1, actual.size());
        assertEquals(TEST_SUPP_DATA1, actual.get(0));
    }

    @Test
    public void testFrom_SingleAlignment_AlternativeDelimiter()
    {
        final List<SupplementaryReadData> actual = SupplementaryReadData.extractAlignments(TEST_SUPP_DATA1.asCsv());

        assertNotNull(actual);
        assertEquals(1, actual.size());
        assertEquals(TEST_SUPP_DATA1, actual.get(0));
    }

    @Test
    public void testFrom_TwoAlignments_WithEndingDelim()
    {
        final List<SupplementaryReadData> actual =
                SupplementaryReadData.extractAlignments(TEST_SUPP_DATA1.asSamTag() + ';' + TEST_SUPP_DATA2.asSamTag() + ';');

        assertNotNull(actual);
        assertEquals(2, actual.size());
        assertEquals(TEST_SUPP_DATA1, actual.get(0));
        assertEquals(TEST_SUPP_DATA2, actual.get(1));
    }

    @Test
    public void testFrom_TwoAlignments_WithoutEndingDelim()
    {
        final List<SupplementaryReadData> actual =
                SupplementaryReadData.extractAlignments(TEST_SUPP_DATA1.asSamTag() + ';' + TEST_SUPP_DATA2.asSamTag());

        assertNotNull(actual);
        assertEquals(2, actual.size());
        assertEquals(TEST_SUPP_DATA1, actual.get(0));
        assertEquals(TEST_SUPP_DATA2, actual.get(1));
    }

    @Test
    public void testFrom_TwoAlignment_EmptyFirst_ReturnsNull()
    {
        assertNull(SupplementaryReadData.extractAlignments(";" + TEST_SUPP_DATA1.asSamTag()));
    }

    @Test
    public void testFrom_TwoAlignment_EmptyLast_ReturnsNull()
    {
        assertNull(SupplementaryReadData.extractAlignments(TEST_SUPP_DATA1.asSamTag() + ";;"));
    }

    @Test
    public void testFirstAlignmentFrom_Emtpy_ReturnsNull()
    {
        assertNull(SupplementaryReadData.extractAlignment(""));
    }

    @Test
    public void testFirstAlignmentFrom_SingleAlignment_WithEndingDelim()
    {
        final SupplementaryReadData actual = SupplementaryReadData.extractAlignment(TEST_SUPP_DATA1.asSamTag() + ';');

        assertNotNull(actual);
        assertEquals(TEST_SUPP_DATA1, actual);
    }

    @Test
    public void testFirstAlignmentFrom_SingleAlignment_WithoutEndingDelim()
    {
        final SupplementaryReadData actual = SupplementaryReadData.extractAlignment(TEST_SUPP_DATA1.asSamTag());

        assertNotNull(actual);
        assertEquals(TEST_SUPP_DATA1, actual);
    }

    @Test
    public void testFirstAlignmentFrom_SingleAlignment_AlternativeDelimiter()
    {
        final SupplementaryReadData actual = SupplementaryReadData.extractAlignment(TEST_SUPP_DATA1.asCsv());

        assertNotNull(actual);
        assertEquals(TEST_SUPP_DATA1, actual);
    }

    @Test
    public void testFirstAlignmentFrom_TwoAlignments_WithEndingDelim()
    {
        final SupplementaryReadData actual =
                SupplementaryReadData.extractAlignment(TEST_SUPP_DATA1.asSamTag() + ';' + TEST_SUPP_DATA2.asSamTag() + ';');

        assertNotNull(actual);
        assertEquals(TEST_SUPP_DATA1, actual);
    }

    @Test
    public void testFirstAlignmentFrom_TwoAlignments_WithoutEndingDelim()
    {
        final SupplementaryReadData actual =
                SupplementaryReadData.extractAlignment(TEST_SUPP_DATA1.asSamTag() + ';' + TEST_SUPP_DATA2.asSamTag());

        assertNotNull(actual);
        assertEquals(TEST_SUPP_DATA1, actual);
    }

    @Test
    public void testFirstAlignmentFrom_TwoAlignment_EmptyFirst_ReturnsNull()
    {
        assertNull(SupplementaryReadData.extractAlignment(";" + TEST_SUPP_DATA1.asSamTag()));
    }

    @Test
    public void testFirstAlignmentFrom_TwoAlignment_EmptyLast()
    {
        final SupplementaryReadData actual = SupplementaryReadData.extractAlignment(TEST_SUPP_DATA1.asSamTag() + ";;");

        assertNotNull(actual);
        assertEquals(TEST_SUPP_DATA1, actual);
    }

    @Test
    public void testAlignmentCount_TwoAlignments_EndingWithDelim()
    {
        final String suppData = TEST_SUPP_DATA1.asSamTag() + ';' + TEST_SUPP_DATA2.asSamTag() + ';';

        assertEquals(2, SupplementaryReadData.alignmentCount(suppData));
    }

    @Test
    public void testAlignmentCount_TwoAlignments_EndingWithoutDelim()
    {
        final String suppData = TEST_SUPP_DATA1.asSamTag() + ';' + TEST_SUPP_DATA2.asSamTag();

        assertEquals(2, SupplementaryReadData.alignmentCount(suppData));
    }
}
