package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.SequenceUtils.buildSequence;
import static com.hartwig.hmftools.panelbuilder.SequenceUtils.isDnaSequenceNormal;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

public class SequenceUtilsTest
{
    private final MockRefGenome mRefGenome;

    public SequenceUtilsTest()
    {
        mRefGenome = new MockRefGenome(true);
        mRefGenome.RefGenomeMap = Map.of(
                "1", "AAAAAAAAAA",
                "2", "GGGGGAAAAA"
        );
        // Compute chromosome lengths based on base sequences.
        mRefGenome.ChromosomeLengths =
                mRefGenome.RefGenomeMap.entrySet().stream().collect(
                        Collectors.toMap(Map.Entry::getKey, e -> e.getValue().length()));
    }

    @Test
    public void testBuildSequenceExactRegion()
    {
        SequenceDefinition def = SequenceDefinition.exactRegion(new ChrBaseRegion("1", 1, 10));
        String actual = buildSequence(mRefGenome, def);
        String expected = "AAAAAAAAAA";
        assertEquals(expected, actual);
    }

    @Test
    public void testBuildSequenceVariant()
    {
        SequenceDefinition def = SequenceDefinition.structuralVariant(
                new ChrBaseRegion("1", 1, 10),
                Orientation.FORWARD,
                "GCGCGCGCGC",
                new ChrBaseRegion("2", 1, 10),
                Orientation.REVERSE);
        String actual = buildSequence(mRefGenome, def);
        String expected = "AAAAAAAAAAGCGCGCGCGCTTTTTCCCCC";
        assertEquals(expected, actual);
    }

    @Test
    public void testBuildSequenceOutOfBounds()
    {
        SequenceDefinition def = SequenceDefinition.exactRegion(new ChrBaseRegion("1", 1, 100));
        assertThrows(IllegalArgumentException.class, () -> buildSequence(mRefGenome, def));
    }

    @Test
    public void testIsDnaSequenceNormal()
    {
        assertTrue(isDnaSequenceNormal("AAAAAA"));
        assertTrue(isDnaSequenceNormal("ACGTTGCA"));
        assertFalse(isDnaSequenceNormal("ACGTNTGCA"));
        assertFalse(isDnaSequenceNormal("NAAAAA"));
        assertFalse(isDnaSequenceNormal("TTTTN"));
        assertFalse(isDnaSequenceNormal("XYZ"));
        assertFalse(isDnaSequenceNormal("NNNNN"));
    }
}
