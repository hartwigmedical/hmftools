package com.hartwig.hmftools.serve.extraction.copynumber;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.extraction.catalog.DealWithDriverInconsistentModeAnnotation;
import com.hartwig.hmftools.serve.extraction.util.GeneChecker;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CopyNumberExtractorTest {

    @Test
    public void canExtractCopyNumbersAmp() {
        CopyNumberExtractor copyNumberExtractor = createTestExtractor();
        KnownCopyNumber amp = copyNumberExtractor.extract("AKT1", EventType.AMPLIFICATION, DealWithDriverInconsistentModeAnnotation.IGNORE);

        assertEquals("AKT1", amp.gene());
        assertEquals(CopyNumberType.AMPLIFICATION, amp.type());
    }

    @Test
    public void canFilterAmpOnUnknownGene() {
        CopyNumberExtractor copyNumberExtractor = createTestExtractor();
        assertNull(copyNumberExtractor.extract("NOT-A-GENE", EventType.AMPLIFICATION, DealWithDriverInconsistentModeAnnotation.IGNORE));
    }

    @Test
    public void canExtractCopyNumbersDel() {
        CopyNumberExtractor copyNumberExtractor = createTestExtractor();
        KnownCopyNumber del = copyNumberExtractor.extract("PTEN", EventType.DELETION, DealWithDriverInconsistentModeAnnotation.IGNORE);

        assertEquals("PTEN", del.gene());
        assertEquals(CopyNumberType.DELETION, del.type());
    }

    @Test
    public void canFilterDelOnUnknownGene() {
        CopyNumberExtractor copyNumberExtractor = createTestExtractor();
        assertNull(copyNumberExtractor.extract("NOT-A-GENE", EventType.DELETION, DealWithDriverInconsistentModeAnnotation.IGNORE));
    }

    @NotNull
    private static CopyNumberExtractor createTestExtractor() {
        return new CopyNumberExtractor(new GeneChecker(Sets.newHashSet("PTEN", "AKT1")), Lists.newArrayList(), true);
    }
}