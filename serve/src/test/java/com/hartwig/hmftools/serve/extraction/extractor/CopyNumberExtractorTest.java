package com.hartwig.hmftools.serve.extraction.extractor;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.copynumber.CopyNumberType;
import com.hartwig.hmftools.serve.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.GeneChecker;
import com.hartwig.hmftools.serve.extraction.GeneCheckerTestFactory;

import org.junit.Test;

public class CopyNumberExtractorTest {

    private static final GeneChecker HG19_GENE_CHECKER = GeneCheckerTestFactory.buildForHG19();

    @Test
    public void canExtractCopyNumbersAmp() {
        CopyNumberExtractor copyNumberExtractor = new CopyNumberExtractor(HG19_GENE_CHECKER);
        KnownCopyNumber amp = copyNumberExtractor.extract("AKT1", EventType.AMPLIFICATION);

        assertEquals("AKT1", amp.gene());
        assertEquals(CopyNumberType.AMPLIFICATION, amp.type());
    }

    @Test
    public void canFilterAmpOnUnknownGene() {
        CopyNumberExtractor copyNumberExtractor = new CopyNumberExtractor(HG19_GENE_CHECKER);
        assertNull(copyNumberExtractor.extract("NOT-A-GENE", EventType.AMPLIFICATION));
    }

    @Test
    public void canExtractCopyNumbersDel() {
        CopyNumberExtractor copyNumberExtractor = new CopyNumberExtractor(HG19_GENE_CHECKER);
        KnownCopyNumber del = copyNumberExtractor.extract("PTEN", EventType.DELETION);

        assertEquals("PTEN", del.gene());
        assertEquals(CopyNumberType.DELETION, del.type());
    }

    @Test
    public void canFilterDelOnUnknownGene() {
        CopyNumberExtractor copyNumberExtractor = new CopyNumberExtractor(HG19_GENE_CHECKER);
        assertNull(copyNumberExtractor.extract("NOT-A-GENE", EventType.DELETION));
    }
}