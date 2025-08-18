package com.hartwig.hmftools.panelbuilder.samplevariants;

import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantProbeBuilder.buildMutationProbe;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantProbeBuilder.buildSglProbe;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantProbeBuilder.buildSvProbe;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.panelbuilder.ProbeTarget;

import org.junit.Test;

public class VariantProbeBuilderTest
{
    private static final String CHR = "1";
    private static final String CHR2 = "2";

    @Test
    public void testBuildMutationProbeSnv()
    {
        ProbeTarget actual = buildMutationProbe(CHR, 101, "A", "C", 11);
        ProbeTarget expected = new ProbeTarget(
                new ChrBaseRegion(CHR, 96, 100),
                Orientation.FORWARD,
                "C",
                new ChrBaseRegion(CHR, 102, 106),
                Orientation.FORWARD);
        assertEquals(expected, actual);
    }

    @Test
    public void testBuildMutationProbeDel()
    {
        ProbeTarget actual = buildMutationProbe(CHR, 101, "A", "", 11);
        ProbeTarget expected = new ProbeTarget(
                new ChrBaseRegion(CHR, 96, 100),
                Orientation.FORWARD,
                "",
                new ChrBaseRegion(CHR, 102, 107),
                Orientation.FORWARD);
        assertEquals(expected, actual);
    }

    @Test
    public void testBuildMutationProbeIndel()
    {
        ProbeTarget actual = buildMutationProbe(CHR, 101, "A", "CG", 11);
        ProbeTarget expected = new ProbeTarget(
                new ChrBaseRegion(CHR, 97, 100),
                Orientation.FORWARD,
                "CG",
                new ChrBaseRegion(CHR, 102, 106),
                Orientation.FORWARD);
        assertEquals(expected, actual);
    }

    @Test
    public void testBuildSglProbe()
    {
        ProbeTarget actual = buildSglProbe(CHR, 101, ORIENT_FWD, "CGT", 11);
        ProbeTarget expected = new ProbeTarget(
                new ChrBaseRegion(CHR, 94, 101),
                Orientation.FORWARD,
                "CGT",
                null,
                null);
        assertEquals(expected, actual);
    }

    @Test
    public void testBuildSvProbe()
    {
        ProbeTarget actual = buildSvProbe(
                CHR, 101, ORIENT_FWD,
                CHR2, 201, ORIENT_REV,
                "GGT", 11);
        ProbeTarget expected = new ProbeTarget(
                new ChrBaseRegion(CHR, 98, 101),
                Orientation.FORWARD,
                "GGT",
                new ChrBaseRegion(CHR2, 201, 204),
                Orientation.FORWARD);
        assertEquals(expected, actual);
    }
}
