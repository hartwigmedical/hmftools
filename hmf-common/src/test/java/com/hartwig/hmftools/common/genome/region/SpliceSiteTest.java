package com.hartwig.hmftools.common.genome.region;

import static com.hartwig.hmftools.common.genome.region.HmfTranscriptRegionUtils.isAcceptorPlusThree;
import static com.hartwig.hmftools.common.genome.region.HmfTranscriptRegionUtils.isDonorMinusOne;
import static com.hartwig.hmftools.common.genome.region.HmfTranscriptRegionUtils.isDonorPlusFive;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;

import org.junit.Test;

public class SpliceSiteTest {

    @Test
    public void testDonorMinusOneForward() {
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
    public void testDonorPlusFiveForward() {
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
    public void testAcceptorPlusThreeForward() {
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
    public void testDonorMinusOneReverse() {
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
    public void testDonorPlusFiveReverse() {
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
    public void testAcceptorPlusThreeReverse() {
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
