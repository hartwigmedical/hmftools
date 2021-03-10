package com.hartwig.hmftools.common.genome.region;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;

import org.junit.Test;

public class SpliceSiteTest {

    @Test
    public void testDonorMinusOneForward() {
        final HmfTranscriptRegion transcript = HmfGenePanelSupplier.allGenesMap37().get("AR");
        assertTrue(transcript.isDonorMinusOne(66766604));
        assertTrue(transcript.isDonorMinusOne(66863249));
        assertTrue(transcript.isDonorMinusOne(66905968));
        assertTrue(transcript.isDonorMinusOne(66931531));
        assertTrue(transcript.isDonorMinusOne(66937464));
        assertTrue(transcript.isDonorMinusOne(66941805));
        assertTrue(transcript.isDonorMinusOne(66942826));
        assertFalse(transcript.isDonorMinusOne(66950461));
    }

    @Test
    public void testDonorPlusFiveForward() {
        final HmfTranscriptRegion transcript = HmfGenePanelSupplier.allGenesMap37().get("AR");
        assertTrue(transcript.isDonorPlusFive(66766604 + 5));
        assertTrue(transcript.isDonorPlusFive(66863249 + 5));
        assertTrue(transcript.isDonorPlusFive(66905968 + 5));
        assertTrue(transcript.isDonorPlusFive(66931531 + 5));
        assertTrue(transcript.isDonorPlusFive(66937464 + 5));
        assertTrue(transcript.isDonorPlusFive(66941805 + 5));
        assertTrue(transcript.isDonorPlusFive(66942826 + 5));
        assertFalse(transcript.isDonorPlusFive(66950461 + 5));
    }

    @Test
    public void testAcceptorPlusThreeForward() {
        final HmfTranscriptRegion transcript = HmfGenePanelSupplier.allGenesMap37().get("AR");
        assertFalse(transcript.isAcceptorPlusThree(66764465 - 3));
        assertTrue(transcript.isAcceptorPlusThree(66863098 - 3));
        assertTrue(transcript.isAcceptorPlusThree(66905852 - 3));
        assertTrue(transcript.isAcceptorPlusThree(66931244 - 3));
        assertTrue(transcript.isAcceptorPlusThree(66937320 - 3));
        assertTrue(transcript.isAcceptorPlusThree(66941675 - 3));
        assertTrue(transcript.isAcceptorPlusThree(66942669 - 3));
        assertTrue(transcript.isAcceptorPlusThree(66943528 - 3));
    }

    @Test
    public void testDonorMinusOneReverse() {
        final HmfTranscriptRegion transcript = HmfGenePanelSupplier.allGenesMap37().get("TP53");

        assertFalse(transcript.isDonorMinusOne(7571720));

        assertTrue(transcript.isDonorMinusOne(7573927));
        assertTrue(transcript.isDonorMinusOne(7576853));
        assertTrue(transcript.isDonorMinusOne(7577019));
        assertTrue(transcript.isDonorMinusOne(7577499));
        assertTrue(transcript.isDonorMinusOne(7578177));
        assertTrue(transcript.isDonorMinusOne(7578371));
        assertTrue(transcript.isDonorMinusOne(7579312));
        assertTrue(transcript.isDonorMinusOne(7579700));
        assertTrue(transcript.isDonorMinusOne(7579839));

        assertFalse(transcript.isDonorMinusOne(7590695));
    }

    @Test
    public void testDonorPlusFiveReverse() {
        final HmfTranscriptRegion transcript = HmfGenePanelSupplier.allGenesMap37().get("TP53");

        assertFalse(transcript.isDonorPlusFive(7571715));

        assertTrue(transcript.isDonorPlusFive(7573922));
        assertTrue(transcript.isDonorPlusFive(7576848));
        assertTrue(transcript.isDonorPlusFive(7577014));
        assertTrue(transcript.isDonorPlusFive(7577494));
        assertTrue(transcript.isDonorPlusFive(7578172));
        assertTrue(transcript.isDonorPlusFive(7578366));
        assertTrue(transcript.isDonorPlusFive(7579307));
        assertTrue(transcript.isDonorPlusFive(7579695));
        assertTrue(transcript.isDonorPlusFive(7579834));

        assertFalse(transcript.isDonorPlusFive(7590690));
    }

    @Test
    public void testAcceptorPlusThreeReverse() {
        final HmfTranscriptRegion transcript = HmfGenePanelSupplier.allGenesMap37().get("TP53");

        assertTrue(transcript.isAcceptorPlusThree(7573008 + 3));
        assertTrue(transcript.isAcceptorPlusThree(7574033 + 3));
        assertTrue(transcript.isAcceptorPlusThree(7576926 + 3));
        assertTrue(transcript.isAcceptorPlusThree(7577155 + 3));
        assertTrue(transcript.isAcceptorPlusThree(7577608 + 3));
        assertTrue(transcript.isAcceptorPlusThree(7578289 + 3));
        assertTrue(transcript.isAcceptorPlusThree(7578554 + 3));
        assertTrue(transcript.isAcceptorPlusThree(7579590 + 3));
        assertTrue(transcript.isAcceptorPlusThree(7579721 + 3));

        assertFalse(transcript.isAcceptorPlusThree(7579940 + 3));
        assertFalse(transcript.isAcceptorPlusThree(7590856 + 3));

    }

}
