package com.hartwig.hmftools.common.dnds;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.Map;

import org.junit.Test;

public class RefCdsCvSupplierTest {

    @Test
    public void testReadFromResource() throws IOException {
        final Map<String, RefCdsCv> refCdsCv = RefCdsCvSupplier.refCdsCv();
        assertEquals(14, refCdsCv.get("TP53").synonymousN());
        assertEquals(745, refCdsCv.get("TP53").missenseN());
        assertEquals(168, refCdsCv.get("TP53").nonsenseN());
        assertEquals(70, refCdsCv.get("TP53").spliceN());
        assertEquals(215, refCdsCv.get("TP53").indelN());
    }

}
