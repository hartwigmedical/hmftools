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

    private ProbeFactory mFactory;

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
        ProbeTarget target = ProbeTarget.exactRegion(new ChrBaseRegion("1", 1, 10));
        Optional<Probe> actual = mFactory.createProbe(target, METADATA);
        Optional<Probe> expected = Optional.of(new Probe(target, "AAAAAAAAAA", METADATA, null, null, DEFAULT_PROBE_QUALITY, 0));
        assertEquals(expected, actual);
    }

    @Test
    public void testCreateProbeVariant()
    {
        ProbeTarget target = ProbeTarget.structuralVariant(
                new ChrBaseRegion("1", 1, 10),
                Orientation.FORWARD,
                "GCGCGCGCGC",
                new ChrBaseRegion("2", 1, 10),
                Orientation.REVERSE);
        Optional<Probe> actual = mFactory.createProbe(target, METADATA);
        Optional<Probe> expected =
                Optional.of(new Probe(target, "AAAAAAAAAAGCGCGCGCGCTTTTTCCCCC", METADATA, null, null, DEFAULT_PROBE_QUALITY, 0.5));
        assertEquals(expected, actual);
    }

    @Test
    public void testCreateProbeOutOfBounds()
    {
        ProbeTarget target = ProbeTarget.exactRegion(new ChrBaseRegion("1", 1, 100));
        Optional<Probe> actual = mFactory.createProbe(target, METADATA);
        assertEquals(actual, Optional.empty());
    }
}
