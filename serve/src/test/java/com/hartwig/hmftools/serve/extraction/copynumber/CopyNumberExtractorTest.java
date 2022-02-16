package com.hartwig.hmftools.serve.extraction.copynumber;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.extraction.catalog.DealWithDriverInconsistentModeAnnotation;
import com.hartwig.hmftools.serve.extraction.codon.CodonExtractor;
import com.hartwig.hmftools.serve.extraction.util.GeneChecker;

import org.jetbrains.annotations.NotNull;
import org.junit.Ignore;
import org.junit.Test;

public class CopyNumberExtractorTest {

    @Test
    @Ignore
    public void canCheckFiltering() {
//        CopyNumberExtractor copyNumberExtractor1 = createTestExtractor(DealWithDriverInconsistentModeAnnotation.IGNORE);
//        assertEquals(CopyNumberType.AMPLIFICATION, copyNumberExtractor1.extract("AKT1", EventType.AMPLIFICATION).type());
//
//        CopyNumberExtractor copyNumberExtractor2 = createTestExtractor(DealWithDriverInconsistentModeAnnotation.IGNORE);
//        assertEquals(CopyNumberType.DELETION, copyNumberExtractor2.extract("AKT1", EventType.DELETION).type());
//
//        CopyNumberExtractor copyNumberExtractor3 = createTestExtractor(DealWithDriverInconsistentModeAnnotation.IGNORE);
//        assertEquals(CopyNumberType.AMPLIFICATION, copyNumberExtractor3.extract("AKT1", EventType.DELETION).type());

//        CopyNumberExtractor copyNumberExtractor2 = createTestExtractor(DealWithDriverInconsistentModeAnnotation.WARN_ONLY);
//        copyNumberExtractor2.extract("AKT1", EventType.AMPLIFICATION);
//
//        CopyNumberExtractor copyNumberExtractor3 = createTestExtractor(DealWithDriverInconsistentModeAnnotation.FILTER);
//        copyNumberExtractor3.extract("AKT1", EventType.AMPLIFICATION);

    }

    @Test
    public void canExtractCopyNumbersAmp() {
        CopyNumberExtractor copyNumberExtractor = createTestExtractor(DealWithDriverInconsistentModeAnnotation.IGNORE);
        KnownCopyNumber amp = copyNumberExtractor.extract("AKT1", EventType.AMPLIFICATION);

        assertEquals("AKT1", amp.gene());
        assertEquals(CopyNumberType.AMPLIFICATION, amp.type());
    }

    @Test
    public void canFilterAmpOnUnknownGene() {
        CopyNumberExtractor copyNumberExtractor = createTestExtractor(DealWithDriverInconsistentModeAnnotation.IGNORE);
        assertNull(copyNumberExtractor.extract("NOT-A-GENE", EventType.AMPLIFICATION));
    }

    @Test
    public void canExtractCopyNumbersDel() {
        CopyNumberExtractor copyNumberExtractor = createTestExtractor(DealWithDriverInconsistentModeAnnotation.IGNORE);
        KnownCopyNumber del = copyNumberExtractor.extract("PTEN", EventType.DELETION);

        assertEquals("PTEN", del.gene());
        assertEquals(CopyNumberType.DELETION, del.type());
    }

    @Test
    public void canFilterDelOnUnknownGene() {
        CopyNumberExtractor copyNumberExtractor = createTestExtractor(DealWithDriverInconsistentModeAnnotation.IGNORE);
        assertNull(copyNumberExtractor.extract("NOT-A-GENE", EventType.DELETION));
    }

    @NotNull
    private static CopyNumberExtractor createTestExtractor(@NotNull DealWithDriverInconsistentModeAnnotation annotation) {
        return new CopyNumberExtractor(new GeneChecker(Sets.newHashSet("PTEN", "AKT1")),
                Lists.newArrayList(), annotation);
    }
}