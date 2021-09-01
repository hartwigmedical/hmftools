package com.hartwig.hmftools.orange;

import static org.junit.Assert.assertTrue;

import com.google.common.io.Resources;

import org.junit.Test;

public class OrangeApiTest {

    private static final String TEST_ORANGE_JSON = Resources.getResource("orange/tumor_sample.orange.json").getPath();

    @Test
    public void canReadOrangeReport() {
        // TODO
        assertTrue(true);
//        OrangeReport report = OrangeApi.read(TEST_ORANGE_JSON);
    }
}