package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.DEFAULT_PROBE_QUALITY;

import static org.junit.Assert.assertEquals;

import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

public class ProbeFactoryTest
{
    private static final TargetMetadata METADATA = new TargetMetadata(TargetMetadata.Type.CUSTOM, "test");

    private final ProbeFactory mFactory;

    public ProbeFactoryTest()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap = Map.of(
                "1", "AAAAAAAAAA",
                "2", "GGGGGAAAAA"
        );
        // Compute chromosome lengths based on base sequences.
        refGenome.ChromosomeLengths =
                refGenome.RefGenomeMap.entrySet().stream().collect(
                        Collectors.toMap(Map.Entry::getKey, e -> e.getValue().length()));
        mFactory = new ProbeFactory(refGenome, null, null);
    }

    @Test
    public void testCreateProbeExactRegion()
    {
        SequenceDefinition def = SequenceDefinition.exactRegion(new ChrBaseRegion("1", 1, 10));
        Optional<Probe> actual = mFactory.createProbe(def, METADATA);
        Optional<Probe> expected = Optional.of(new Probe(def, "AAAAAAAAAA", METADATA, null, null, DEFAULT_PROBE_QUALITY, 0));
        assertEquals(expected, actual);
    }

    @Test
    public void testCreateProbeVariant()
    {
        SequenceDefinition def = SequenceDefinition.structuralVariant(
                new ChrBaseRegion("1", 1, 10),
                Orientation.FORWARD,
                "GCGCGCGCGC",
                new ChrBaseRegion("2", 1, 10),
                Orientation.REVERSE);
        Optional<Probe> actual = mFactory.createProbe(def, METADATA);
        Optional<Probe> expected =
                Optional.of(new Probe(def, "AAAAAAAAAAGCGCGCGCGCTTTTTCCCCC", METADATA, null, null, DEFAULT_PROBE_QUALITY, 0.5));
        assertEquals(expected, actual);
    }

    @Test
    public void testCreateProbeOutOfBounds()
    {
        SequenceDefinition def = SequenceDefinition.exactRegion(new ChrBaseRegion("1", 1, 100));
        Optional<Probe> actual = mFactory.createProbe(def, METADATA);
        assertEquals(actual, Optional.empty());
    }
}
