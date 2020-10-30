package com.hartwig.hmftools.common.pharmacogenetics;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class PGXCallsFileTest {

    private static final String FILE = Resources.getResource("pharmacogenetics/sample_calls.txt").getPath();

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
    public void loadPgxCallsFile() throws IOException {
        List<PGXCalls> pgxCalls = PGXCallsFile.read(FILE);

        assertEquals(18, pgxCalls.size());

        PGXCalls pgxCalls1 = pgxCalls.get(0);
        assertEquals(GENE_1, pgxCalls1.gene());
        assertEquals(POSITION_GRCH37_1, pgxCalls1.positionGRCh37());
        assertEquals(REF_GRCH37_1, pgxCalls1.refGRCh37());
        assertEquals(ALT_GRCH37_1, pgxCalls1.altGRCh37());
        assertEquals(POSITION_GRCH38_1, pgxCalls1.positionGRCh38());
        assertEquals(REF_GRCH38_1, pgxCalls1.refGRCh38());
        assertEquals(ALT_GRCH38_1, pgxCalls1.altGRCh38());
        assertEquals(RSID_1, pgxCalls1.rsid());
        assertEquals(VARIANT_ANNOTATION_1, pgxCalls1.variantAnnotation());
        assertEquals(FILTER_1, pgxCalls1.filter());

        PGXCalls pgxCalls2 = pgxCalls.get(15);
        assertEquals(GENE_2, pgxCalls2.gene());
        assertEquals(POSITION_GRCH37_2, pgxCalls2.positionGRCh37());
        assertEquals(REF_GRCH37_2, pgxCalls2.refGRCh37());
        assertEquals(ALT_GRCH37_2, pgxCalls2.altGRCh37());
        assertEquals(POSITION_GRCH38_2, pgxCalls2.positionGRCh38());
        assertEquals(REF_GRCH38_2, pgxCalls2.refGRCh38());
        assertEquals(ALT_GRCH38_2, pgxCalls2.altGRCh38());
        assertEquals(RSID_2, pgxCalls2.rsid());
        assertEquals(VARIANT_ANNOTATION_2, pgxCalls2.variantAnnotation());
        assertEquals(FILTER_2, pgxCalls2.filter());
    }

}