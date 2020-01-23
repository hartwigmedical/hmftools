package com.hartwig.hmftools.protect.actionability.cancertype;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import org.junit.Test;

public class PrimaryTumorToDOIDMapperTest {

    @Test
    public void canCreateProductionResource() throws IOException {
        assertNotNull(PrimaryTumorToDOIDMapper.createFromResource());
    }
}
