package com.hartwig.hmftools.orange.algo.sigs;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.Map;

import com.google.common.io.Resources;

import org.junit.Test;

public class SigsEtiologiesLoaderTest
{
    private static final String SIGNATURES_ETIOLOGY_TSV =
            Resources.getResource("test_run_resources/sigs/signatures_etiology.tsv").getPath();

    @Test
    public void canReadSignaturesEtiologyTsv() throws IOException
    {
        Map<String, String> etiologyPerSignature = SigsEtiologiesLoader.read(SIGNATURES_ETIOLOGY_TSV);
        assertEquals(30, etiologyPerSignature.size());
    }
}