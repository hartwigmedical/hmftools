package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.genome.region.Strand.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Strand.REVERSE;
import static com.hartwig.hmftools.common.variant.SpliceSites.getAcceptorPosition;
import static com.hartwig.hmftools.common.variant.SpliceSites.getDonorPosition;
import static com.hartwig.hmftools.common.variant.SpliceSites.isAcceptorPlusThree;
import static com.hartwig.hmftools.common.variant.SpliceSites.isDonorMinusOne;
import static com.hartwig.hmftools.common.variant.SpliceSites.isDonorPlusFive;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;

import org.junit.Test;

public class SpliceSiteTest
{
    @Test
    public void testDonorSites()
    {
        assertEquals(-2, getDonorPosition(99, 100, FORWARD));
        assertEquals(-1, getDonorPosition(100, 100, FORWARD));
        assertEquals(1, getDonorPosition(101, 100, FORWARD));
        assertEquals(2, getDonorPosition(102, 100, FORWARD));
        assertEquals(5, getDonorPosition(105, 100, FORWARD));

        assertEquals(-2, getDonorPosition(101, 100, REVERSE));
        assertEquals(-1, getDonorPosition(100, 100, REVERSE));
        assertEquals(1, getDonorPosition(99, 100, REVERSE));
        assertEquals(2, getDonorPosition(98, 100, REVERSE));
        assertEquals(5, getDonorPosition(95, 100, REVERSE));
    }

    @Test
    public void testAcceptorSites()
    {
        assertEquals(3, getAcceptorPosition(97, 100, FORWARD));
        assertEquals(2, getAcceptorPosition(98, 100, FORWARD));
        assertEquals(1, getAcceptorPosition(99, 100, FORWARD));
        assertEquals(-1, getAcceptorPosition(100, 100, FORWARD));

        assertEquals(3, getAcceptorPosition(103, 100, REVERSE));
        assertEquals(2, getAcceptorPosition(102, 100, REVERSE));
        assertEquals(1, getAcceptorPosition(101, 100, REVERSE));
        assertEquals(-1, getAcceptorPosition(100, 100, REVERSE));
    }

    // TODO - simplify, so unnecessary
    @Test
    public void testDonorMinusOneForward()
    {
        final HmfTranscriptRegion transcript = HmfGenePanelSupplier.allGenesMap37().get("AR");
        assertTrue(isDonorMinusOne(transcript, 66766604));
        assertTrue(isDonorMinusOne(transcript, 66863249));
        assertTrue(isDonorMinusOne(transcript, 66905968));
        assertTrue(isDonorMinusOne(transcript, 66931531));
        assertTrue(isDonorMinusOne(transcript, 66937464));
        assertTrue(isDonorMinusOne(transcript, 66941805));
        assertTrue(isDonorMinusOne(transcript, 66942826));
        assertFalse(isDonorMinusOne(transcript, 66950461));
    }

    @Test
    public void testDonorPlusFiveForward()
    {
        final HmfTranscriptRegion transcript = HmfGenePanelSupplier.allGenesMap37().get("AR");
        assertTrue(isDonorPlusFive(transcript, 66766604 + 5));
        assertTrue(isDonorPlusFive(transcript, 66863249 + 5));
        assertTrue(isDonorPlusFive(transcript, 66905968 + 5));
        assertTrue(isDonorPlusFive(transcript, 66931531 + 5));
        assertTrue(isDonorPlusFive(transcript, 66937464 + 5));
        assertTrue(isDonorPlusFive(transcript, 66941805 + 5));
        assertTrue(isDonorPlusFive(transcript, 66942826 + 5));
        assertFalse(isDonorPlusFive(transcript, 66950461 + 5));
    }

    @Test
    public void testAcceptorPlusThreeForward()
    {
        final HmfTranscriptRegion transcript = HmfGenePanelSupplier.allGenesMap37().get("AR");
        assertFalse(isAcceptorPlusThree(transcript, 66764465 - 3));
        assertTrue(isAcceptorPlusThree(transcript, 66863098 - 3));
        assertTrue(isAcceptorPlusThree(transcript, 66905852 - 3));
        assertTrue(isAcceptorPlusThree(transcript, 66931244 - 3));
        assertTrue(isAcceptorPlusThree(transcript, 66937320 - 3));
        assertTrue(isAcceptorPlusThree(transcript, 66941675 - 3));
        assertTrue(isAcceptorPlusThree(transcript, 66942669 - 3));
        assertTrue(isAcceptorPlusThree(transcript, 66943528 - 3));
    }

    @Test
    public void testDonorMinusOneReverse()
    {
        final HmfTranscriptRegion transcript = HmfGenePanelSupplier.allGenesMap37().get("TP53");

        assertFalse(isDonorMinusOne(transcript, 7571720));

        assertTrue(isDonorMinusOne(transcript, 7573927));
        assertTrue(isDonorMinusOne(transcript, 7576853));
        assertTrue(isDonorMinusOne(transcript, 7577019));
        assertTrue(isDonorMinusOne(transcript, 7577499));
        assertTrue(isDonorMinusOne(transcript, 7578177));
        assertTrue(isDonorMinusOne(transcript, 7578371));
        assertTrue(isDonorMinusOne(transcript, 7579312));
        assertTrue(isDonorMinusOne(transcript, 7579700));
        assertTrue(isDonorMinusOne(transcript, 7579839));

        assertFalse(isDonorMinusOne(transcript, 7590695));
    }

    @Test
    public void testDonorPlusFiveReverse()
    {
        final HmfTranscriptRegion transcript = HmfGenePanelSupplier.allGenesMap37().get("TP53");

        assertFalse(isDonorPlusFive(transcript, 7571715));

        assertTrue(isDonorPlusFive(transcript, 7573922));
        assertTrue(isDonorPlusFive(transcript, 7576848));
        assertTrue(isDonorPlusFive(transcript, 7577014));
        assertTrue(isDonorPlusFive(transcript, 7577494));
        assertTrue(isDonorPlusFive(transcript, 7578172));
        assertTrue(isDonorPlusFive(transcript, 7578366));
        assertTrue(isDonorPlusFive(transcript, 7579307));
        assertTrue(isDonorPlusFive(transcript, 7579695));
        assertTrue(isDonorPlusFive(transcript, 7579834));

        assertFalse(isDonorPlusFive(transcript, 7590690));
    }

    @Test
    public void testAcceptorPlusThreeReverse()
    {
        final HmfTranscriptRegion transcript = HmfGenePanelSupplier.allGenesMap37().get("TP53");

        assertTrue(isAcceptorPlusThree(transcript, 7573008 + 3));
        assertTrue(isAcceptorPlusThree(transcript, 7574033 + 3));
        assertTrue(isAcceptorPlusThree(transcript, 7576926 + 3));
        assertTrue(isAcceptorPlusThree(transcript, 7577155 + 3));
        assertTrue(isAcceptorPlusThree(transcript, 7577608 + 3));
        assertTrue(isAcceptorPlusThree(transcript, 7578289 + 3));
        assertTrue(isAcceptorPlusThree(transcript, 7578554 + 3));
        assertTrue(isAcceptorPlusThree(transcript, 7579590 + 3));
        assertTrue(isAcceptorPlusThree(transcript, 7579721 + 3));

        assertFalse(isAcceptorPlusThree(transcript, 7579940 + 3));
        assertFalse(isAcceptorPlusThree(transcript, 7590856 + 3));

    }

}
