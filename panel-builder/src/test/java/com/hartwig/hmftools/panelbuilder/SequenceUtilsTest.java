package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.SequenceUtils.buildIndelProbe;
import static com.hartwig.hmftools.panelbuilder.SequenceUtils.buildSequence;
import static com.hartwig.hmftools.panelbuilder.SequenceUtils.buildSglProbe;
import static com.hartwig.hmftools.panelbuilder.SequenceUtils.buildSvProbe;
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
    public void testBuildSequenceSingleRegion()
    {
        SequenceDefinition def = SequenceDefinition.singleRegion(new ChrBaseRegion("1", 1, 10));
        String actual = buildSequence(mRefGenome, def);
        String expected = "AAAAAAAAAA";
        assertEquals(expected, actual);
    }

    @Test
    public void testBuildSequenceVariant()
    {
        SequenceDefinition def = new SequenceDefinition(
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
        SequenceDefinition def = SequenceDefinition.singleRegion(new ChrBaseRegion("1", 1, 100));
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

    @Test
    public void testBuildIndelProbeSnv()
    {
        SequenceDefinition actual = buildIndelProbe("1", 101, "A", "C", 11);
        SequenceDefinition expected = new SequenceDefinition(
                new ChrBaseRegion("1", 96, 100),
                Orientation.FORWARD,
                "C",
                new ChrBaseRegion("1", 102, 106),
                Orientation.FORWARD);
        assertEquals(expected, actual);
    }

    @Test
    public void testBuildIndelProbeDel()
    {
        SequenceDefinition actual = buildIndelProbe("1", 101, "A", "", 11);
        SequenceDefinition expected = new SequenceDefinition(
                new ChrBaseRegion("1", 96, 100),
                Orientation.FORWARD,
                "",
                new ChrBaseRegion("1", 102, 107),
                Orientation.FORWARD);
        assertEquals(expected, actual);
    }

    @Test
    public void testBuildIndelProbeIndel()
    {
        SequenceDefinition actual = buildIndelProbe("1", 101, "A", "CG", 11);
        SequenceDefinition expected = new SequenceDefinition(
                new ChrBaseRegion("1", 97, 100),
                Orientation.FORWARD,
                "CG",
                new ChrBaseRegion("1", 102, 106),
                Orientation.FORWARD);
        assertEquals(expected, actual);
    }

    @Test
    public void testBuildSglProbe()
    {
        SequenceDefinition actual = buildSglProbe("1", 101, Orientation.FORWARD, "CGT", 11);
        SequenceDefinition expected = new SequenceDefinition(
                new ChrBaseRegion("1", 94, 101),
                Orientation.FORWARD,
                "CGT",
                null,
                null);
        assertEquals(expected, actual);
    }

    @Test
    public void testBuildSvProbe()
    {
        SequenceDefinition actual = buildSvProbe(
                "1", 101, Orientation.FORWARD,
                "2", 201, Orientation.REVERSE,
                "GGT", 11);
        SequenceDefinition expected = new SequenceDefinition(
                new ChrBaseRegion("1", 98, 101),
                Orientation.FORWARD,
                "GGT",
                new ChrBaseRegion("2", 201, 204),
                Orientation.FORWARD);
        assertEquals(expected, actual);
    }
}
