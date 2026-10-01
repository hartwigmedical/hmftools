package com.hartwig.hmftools.viridian.integration;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.viridian.integration.align.ViralInsertAlignment;
import com.hartwig.hmftools.viridian.integration.extract.BreakendSupport;
import com.hartwig.hmftools.viridian.integration.extract.HostBreakend;
import com.hartwig.hmftools.viridian.integration.extract.HostVariantCandidate;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.ViralContig;

import org.junit.Test;

public class IntegrationTest
{
    private static final ViralContig CONTIG = new ViralContig("v1", 7906, "Virus v1", new OncologyGroup("Group A"));

    // An insert that aligned nowhere is not an integration, however the host variant looks.
    @Test
    public void testNotPlausibleWhenInsertAlignedNowhere()
    {
        Integration integration = new Integration(candidate(), null);
        assertFalse(integration.isAligned());
        assertFalse(integration.isPlausible());
    }

    @Test
    public void testPlausibleWhenAlignmentClearsThresholds()
    {
        Integration integration = new Integration(candidate(), alignment(40));
        assertTrue(integration.isAligned());
        assertTrue(integration.isPlausible());
    }

    @Test
    public void testNotPlausibleWhenAlignmentDoesNotClearThresholds()
    {
        Integration integration = new Integration(candidate(), alignment(20));
        assertTrue(integration.isAligned());
        assertFalse(integration.isPlausible());
    }

    private static HostVariantCandidate candidate()
    {
        HostBreakend breakend = new HostBreakend(
                "sgl_1", new BasePosition("chr1", 1000), FORWARD, new BreakendSupport(12, 60, 0.25, 0, 40));
        return new HostVariantCandidate(
                StructuralVariantType.SGL, "PASS", breakend, null, "ACGTACGTACGTACGTACGTACGTACGTAC", false, null, "");
    }

    private static ViralInsertAlignment alignment(int alignerScore)
    {
        return new ViralInsertAlignment(CONTIG, 1500, FORWARD, cigarFromStr("50M"), alignerScore, 1, 50);
    }
}
