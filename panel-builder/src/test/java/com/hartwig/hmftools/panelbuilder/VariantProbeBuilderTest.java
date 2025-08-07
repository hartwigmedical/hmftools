package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantProbeBuilder.buildMutationProbe;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantProbeBuilder.buildSglProbe;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.panelbuilder.samplevariants.VariantProbeData;

import org.junit.Test;

public class VariantProbeBuilderTest
{
    private static final String CHR = "1";

    private final MockRefGenome mRefGenome;

    public VariantProbeBuilderTest()
    {
        mRefGenome = new MockRefGenome(true);
        mRefGenome.RefGenomeMap = Map.of(CHR, "A".repeat(200));
    }

    @Test
    public void testBuildMutationProbeSnv()
    {
        VariantProbeData actual = buildMutationProbe(CHR, 101, "A", "C", 11, mRefGenome);
        VariantProbeData expected = new VariantProbeData(
                "AAAAACAAAAA",
                new ChrBaseRegion(CHR, 96, 100),
                "C",
                new ChrBaseRegion(CHR, 102, 106));
        assertEquals(expected, actual);
    }

    @Test
    public void testBuildMutationProbeDel()
    {
        VariantProbeData actual = buildMutationProbe(CHR, 101, "A", "", 11, mRefGenome);
        VariantProbeData expected = new VariantProbeData(
                "AAAAAAAAAAA",
                new ChrBaseRegion(CHR, 96, 100),
                "",
                new ChrBaseRegion(CHR, 102, 107));
        assertEquals(expected, actual);
    }

    @Test
    public void testBuildMutationProbeIndel()
    {
        VariantProbeData actual = buildMutationProbe(CHR, 101, "A", "CG", 11, mRefGenome);
        VariantProbeData expected = new VariantProbeData(
                "AAAACGAAAAA",
                new ChrBaseRegion(CHR, 97, 100),
                "CG",
                new ChrBaseRegion(CHR, 102, 106));
        assertEquals(expected, actual);
    }

    @Test
    public void testBuildSglProbe()
    {
        VariantProbeData actual = buildSglProbe(CHR, 101, ORIENT_FWD, "CGT", 11, mRefGenome);
        VariantProbeData expected = new VariantProbeData(
                "AAAAAAAACGT",
                new ChrBaseRegion(CHR, 94, 101),
                "CGT",
                null);
        assertEquals(expected, actual);
    }

    // TODO: more tests
}
