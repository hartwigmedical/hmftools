package com.hartwig.hmftools.serve;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import org.junit.Test;

public class ServeLocalConfigProviderTest {

    @Test
    public void canCreateLocalConfig() throws IOException {
        assertNotNull(ServeLocalConfigProvider.create());
    }
}