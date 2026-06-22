package com.hartwig.hmftools.esvee.common.saga;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.region.BasePosition;

import org.junit.Test;

public class SagaAssemblyTest
{
    private static final String VALID_SEQUENCE =
            "AAATTTGGGGCCCCAAAAATTTTTGGGGGGCCCCCCAAAAAAAATTTTTTTGGGGGGGGCCCCCCCCAAAAAAAAA";

    private static final String MINIMAL_SEQUENCE = "AAAAAA";

    @Test
    public void testValidInsertionWith2JunctionOffsets()
    {
        String label = "TestInsert|chr1:100:1|chr1:101:-1|50|75";

        SagaAssembly assembly = SagaAssembly.fromFastaRecord(label, VALID_SEQUENCE);

        SagaVariant expectedVariant = new SagaVariant(
                "TestInsert",
                new SagaBreakend(new BasePosition("chr1", 100), FORWARD),
                new SagaBreakend(new BasePosition("chr1", 101), REVERSE),
                VALID_SEQUENCE.substring(50, 75)
        );

        assertEquals(expectedVariant, assembly.variant());
        assertEquals(List.of(50, 75), assembly.junctionOffsets());
        assertEquals(VALID_SEQUENCE, assembly.sequence());
        assertTrue(assembly.variant().isInsert());
    }

    @Test
    public void testValidDeletionWith1JunctionOffset()
    {
        String label = "TestDelete|chr1:100:1|chr1:150:-1|50";

        SagaAssembly assembly = SagaAssembly.fromFastaRecord(label, VALID_SEQUENCE);

        SagaVariant expectedVariant = new SagaVariant(
                "TestDelete",
                new SagaBreakend(new BasePosition("chr1", 100), FORWARD),
                new SagaBreakend(new BasePosition("chr1", 150), REVERSE),
                ""
        );

        assertEquals(expectedVariant, assembly.variant());
        assertEquals(List.of(50), assembly.junctionOffsets());
        assertEquals(VALID_SEQUENCE, assembly.sequence());
        assertFalse(assembly.variant().isInsert());
    }

    @Test
    public void testInvalidLabelTooFewParts()
    {
        String invalidLabel = "TestInsert|chr1:100:1|chr1:101:-1";

        assertThrows(IllegalArgumentException.class,
                () -> SagaAssembly.fromFastaRecord(invalidLabel, VALID_SEQUENCE));
    }

    @Test
    public void testInvalidLabelTooManyParts()
    {
        String invalidLabel = "TestInsert|chr1:100:1|chr1:101:-1|50|75|extra";

        assertThrows(IllegalArgumentException.class,
                () -> SagaAssembly.fromFastaRecord(invalidLabel, VALID_SEQUENCE));
    }

    @Test
    public void testJunctionOffsetOutOfBoundsZero()
    {
        String label = "TestInsert|chr1:100:1|chr1:101:-1|0|75";

        assertThrows(IllegalArgumentException.class,
                () -> SagaAssembly.fromFastaRecord(label, VALID_SEQUENCE));
    }

    @Test
    public void testJunctionOffsetOutOfBoundsAtSequenceLength()
    {
        String label = "TestInsert|chr1:100:1|chr1:101:-1|50|" + VALID_SEQUENCE.length();

        assertThrows(IllegalArgumentException.class,
                () -> SagaAssembly.fromFastaRecord(label, VALID_SEQUENCE));
    }

    @Test
    public void testJunctionOffsetOutOfBoundsBeyondSequenceLength()
    {
        String label = "TestInsert|chr1:100:1|chr1:101:-1|50|" + (VALID_SEQUENCE.length() + 10);

        assertThrows(Exception.class,
                () -> SagaAssembly.fromFastaRecord(label, VALID_SEQUENCE));
    }

    @Test
    public void testJunctionOffsetsNotInAscendingOrder()
    {
        String label = "TestInsert|chr1:100:1|chr1:101:-1|75|50";

        assertThrows(Exception.class,
                () -> SagaAssembly.fromFastaRecord(label, VALID_SEQUENCE));
    }

    @Test
    public void testJunctionOffsetsEqual()
    {
        String label = "TestInsert|chr1:100:1|chr1:101:-1|50|50";

        assertThrows(IllegalArgumentException.class,
                () -> SagaAssembly.fromFastaRecord(label, VALID_SEQUENCE));
    }

    @Test
    public void testBreakendDifferentChromosomes()
    {
        String label = "TestInsert|chr1:100:1|chr2:101:-1|50|75";

        assertThrows(IllegalArgumentException.class,
                () -> SagaAssembly.fromFastaRecord(label, VALID_SEQUENCE));
    }

    @Test
    public void testBreakendWrongOrientationBothForward()
    {
        String label = "TestInsert|chr1:100:1|chr1:101:1|50|75";

        assertThrows(IllegalArgumentException.class,
                () -> SagaAssembly.fromFastaRecord(label, VALID_SEQUENCE));
    }

    @Test
    public void testBreakendWrongOrientationBothReverse()
    {
        String label = "TestInsert|chr1:100:-1|chr1:101:-1|50|75";

        assertThrows(IllegalArgumentException.class,
                () -> SagaAssembly.fromFastaRecord(label, VALID_SEQUENCE));
    }

    @Test
    public void testBreakendReversePositions()
    {
        String label = "TestInsert|chr1:101:1|chr1:100:-1|50|75";

        assertThrows(IllegalArgumentException.class,
                () -> SagaAssembly.fromFastaRecord(label, VALID_SEQUENCE));
    }

    @Test
    public void testBreakendEqualPositions()
    {
        String label = "TestInsert|chr1:100:1|chr1:100:-1|50|75";

        assertThrows(IllegalArgumentException.class,
                () -> SagaAssembly.fromFastaRecord(label, VALID_SEQUENCE));
    }

    @Test
    public void testBreakendInvalidFormat()
    {
        String label = "TestInsert|chr1:1|chr1:101:-1|50|75";

        assertThrows(IllegalArgumentException.class,
                () -> SagaAssembly.fromFastaRecord(label, VALID_SEQUENCE));
    }

    @Test
    public void testInsertionMinimalBoundaries()
    {
        String label = "TestMin|chr1:100:1|chr1:101:-1|1|2";

        SagaAssembly assembly = SagaAssembly.fromFastaRecord(label, MINIMAL_SEQUENCE);

        assertEquals(List.of(1, 2), assembly.junctionOffsets());
        assertEquals("A", assembly.variant().insertSequence());
        assertTrue(assembly.variant().isInsert());
    }

    @Test
    public void testInsertionMaximalJunctionOffsets()
    {
        String label = "TestMax|chr1:100:1|chr1:101:-1|1|" + (VALID_SEQUENCE.length() - 1);

        SagaAssembly assembly = SagaAssembly.fromFastaRecord(label, VALID_SEQUENCE);

        assertEquals(List.of(1, VALID_SEQUENCE.length() - 1), assembly.junctionOffsets());
        assertEquals(VALID_SEQUENCE.substring(1, VALID_SEQUENCE.length() - 1),
                assembly.variant().insertSequence());
    }

    @Test
    public void testDeletionMinimalJunctionOffset()
    {
        String label = "TestDelMin|chr1:100:1|chr1:150:-1|1";

        SagaAssembly assembly = SagaAssembly.fromFastaRecord(label, MINIMAL_SEQUENCE);

        assertEquals(List.of(1), assembly.junctionOffsets());
        assertEquals("", assembly.variant().insertSequence());
        assertFalse(assembly.variant().isInsert());
    }

    @Test
    public void testDeletionMaximalJunctionOffset()
    {
        String label = "TestDelMax|chr1:100:1|chr1:150:-1|" + (VALID_SEQUENCE.length() - 1);

        SagaAssembly assembly = SagaAssembly.fromFastaRecord(label, VALID_SEQUENCE);

        assertEquals(List.of(VALID_SEQUENCE.length() - 1), assembly.junctionOffsets());
        assertEquals("", assembly.variant().insertSequence());
    }
}