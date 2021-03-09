package com.hartwig.hmftools.common.genome.region;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;

import org.junit.Test;

public class SpliceSiteTest {

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
    public void testAcceptorPlusThreeReverse() {
        final HmfTranscriptRegion transcript = HmfGenePanelSupplier.allGenesMap37().get("TP53");

        assertTrue(transcript.isAcceptorPlusThree(7573008 + 4));
        assertTrue(transcript.isAcceptorPlusThree(7574033 + 4));
        assertTrue(transcript.isAcceptorPlusThree(7576926 + 4));
        assertTrue(transcript.isAcceptorPlusThree(7577155 + 4));
        assertTrue(transcript.isAcceptorPlusThree(7577608 + 4));
        assertTrue(transcript.isAcceptorPlusThree(7578289 + 4));
        assertTrue(transcript.isAcceptorPlusThree(7578554 + 4));
        assertTrue(transcript.isAcceptorPlusThree(7579590 + 4));
        assertTrue(transcript.isAcceptorPlusThree(7579721 + 4));

        assertFalse(transcript.isAcceptorPlusThree(7579940 + 4));
        assertFalse(transcript.isAcceptorPlusThree(7590856 + 4));

    }

}
