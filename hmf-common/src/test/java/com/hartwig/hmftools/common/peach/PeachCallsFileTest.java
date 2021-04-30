package com.hartwig.hmftools.common.peach;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class PeachCallsFileTest {

    private static final String FILE = Resources.getResource("peach/sample_calls.txt").getPath();

    private static final String GENE_1 = "GENE";
    private static final String POSITION_GRCH37_1 = "1:9740";
    private static final String REF_GRCH37_1 = "A";
    private static final String ALT_GRCH37_1 = "A";
    private static final String POSITION_GRCH38_1 = "1:9750";
    private static final String ALT_GRCH38_1 = "A";
    private static final String REF_GRCH38_1 = "A";
    private static final String RSID_1 = "rs456";
    private static final String VARIANT_ANNOTATION_1 = "REF_CALL";
    private static final String FILTER_1 = "NO_CALL";

    private static final String GENE_2 = "GENE";
    private static final String POSITION_GRCH37_2 = "1:8728";
    private static final String REF_GRCH37_2 = "GCCGT";
    private static final String ALT_GRCH37_2 = "GCCGT";
    private static final String POSITION_GRCH38_2 = "1:8747";
    private static final String ALT_GRCH38_2 = "GCCGT";
    private static final String REF_GRCH38_2 = "GCCGT";
    private static final String RSID_2 = "rs9842";
    private static final String VARIANT_ANNOTATION_2 = "1056G>A";
    private static final String FILTER_2 = "PASS";

    @Test
    public void loadPeachCallsFile() throws IOException {
        List<PeachCalls> peachCalls = PeachCallsFile.read(FILE);

        assertEquals(18, peachCalls.size());

        PeachCalls peachCalls1 = peachCalls.get(0);
        assertEquals(GENE_1, peachCalls1.gene());
        assertEquals(POSITION_GRCH37_1, peachCalls1.positionGRCh37());
        assertEquals(REF_GRCH37_1, peachCalls1.refGRCh37());
        assertEquals(ALT_GRCH37_1, peachCalls1.altGRCh37());
        assertEquals(POSITION_GRCH38_1, peachCalls1.positionGRCh38());
        assertEquals(REF_GRCH38_1, peachCalls1.refGRCh38());
        assertEquals(ALT_GRCH38_1, peachCalls1.altGRCh38());
        assertEquals(RSID_1, peachCalls1.rsid());
        assertEquals(VARIANT_ANNOTATION_1, peachCalls1.variantAnnotation());
        assertEquals(FILTER_1, peachCalls1.filter());

        PeachCalls peachCalls2 = peachCalls.get(15);
        assertEquals(GENE_2, peachCalls2.gene());
        assertEquals(POSITION_GRCH37_2, peachCalls2.positionGRCh37());
        assertEquals(REF_GRCH37_2, peachCalls2.refGRCh37());
        assertEquals(ALT_GRCH37_2, peachCalls2.altGRCh37());
        assertEquals(POSITION_GRCH38_2, peachCalls2.positionGRCh38());
        assertEquals(REF_GRCH38_2, peachCalls2.refGRCh38());
        assertEquals(ALT_GRCH38_2, peachCalls2.altGRCh38());
        assertEquals(RSID_2, peachCalls2.rsid());
        assertEquals(VARIANT_ANNOTATION_2, peachCalls2.variantAnnotation());
        assertEquals(FILTER_2, peachCalls2.filter());
    }

}