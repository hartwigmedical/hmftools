package com.hartwig.hmftools.serve.extraction.copynumber;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.DriverGeneTestFactory;
import com.hartwig.hmftools.serve.extraction.util.DriverInconsistencyMode;
import com.hartwig.hmftools.serve.extraction.util.GeneChecker;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CopyNumberExtractorTest {

    @Test
    public void canCheckFilterInCatalog() {
        CopyNumberExtractor copyNumberExtractorIgnore = createTestExtractor(DriverInconsistencyMode.IGNORE);
        assertEquals(CopyNumberType.AMPLIFICATION, copyNumberExtractorIgnore.extract("AKT1", EventType.AMPLIFICATION).type());

        CopyNumberExtractor copyNumberExtractorWarn = createTestExtractor(DriverInconsistencyMode.WARN_ONLY);
        assertEquals(CopyNumberType.AMPLIFICATION, copyNumberExtractorWarn.extract("AKT1", EventType.AMPLIFICATION).type());

        CopyNumberExtractor copyNumberExtractorFilter = createTestExtractor(DriverInconsistencyMode.FILTER);
        assertEquals(CopyNumberType.OVER_EXPRESSION, copyNumberExtractorFilter.extract("AKT1", EventType.OVER_EXPRESSION).type());

        CopyNumberExtractor copyNumberExtractorFilterDel = createTestExtractor(DriverInconsistencyMode.FILTER);
        assertNull(copyNumberExtractorFilterDel.extract("AKT1", EventType.DELETION));
    }

    @Test
    public void canCheckFilterNotInCatalog() {
        CopyNumberExtractor copyNumberExtractorIgnore = createTestExtractor(DriverInconsistencyMode.IGNORE);
        assertEquals(CopyNumberType.AMPLIFICATION, copyNumberExtractorIgnore.extract("KRAS", EventType.AMPLIFICATION).type());

        CopyNumberExtractor copyNumberExtractorWarn = createTestExtractor(DriverInconsistencyMode.WARN_ONLY);
        assertEquals(CopyNumberType.AMPLIFICATION, copyNumberExtractorWarn.extract("PTEN", EventType.AMPLIFICATION).type());

        CopyNumberExtractor copyNumberExtractorFilter = createTestExtractor(DriverInconsistencyMode.FILTER);
        assertNull(copyNumberExtractorFilter.extract("PTEN", EventType.AMPLIFICATION));
    }

    @Test
    public void canExtractCopyNumbersAmp() {
        CopyNumberExtractor copyNumberExtractor = createTestExtractor(DriverInconsistencyMode.IGNORE);
        KnownCopyNumber amp = copyNumberExtractor.extract("AKT1", EventType.AMPLIFICATION);

        assertEquals("AKT1", amp.gene());
        assertEquals(CopyNumberType.AMPLIFICATION, amp.type());
    }

    @Test
    public void canFilterAmpOnUnknownGene() {
        CopyNumberExtractor copyNumberExtractor = createTestExtractor(DriverInconsistencyMode.IGNORE);
        assertNull(copyNumberExtractor.extract("NOT-A-GENE", EventType.AMPLIFICATION));
    }

    @Test
    public void canExtractCopyNumbersDel() {
        CopyNumberExtractor copyNumberExtractor = createTestExtractor(DriverInconsistencyMode.IGNORE);
        KnownCopyNumber del = copyNumberExtractor.extract("PTEN", EventType.UNDER_EXPRESSION);

        assertEquals("PTEN", del.gene());
        assertEquals(CopyNumberType.UNDER_EXPRESSION, del.type());
    }

    @Test
    public void canFilterDelOnUnknownGene() {
        CopyNumberExtractor copyNumberExtractor = createTestExtractor(DriverInconsistencyMode.IGNORE);
        assertNull(copyNumberExtractor.extract("NOT-A-GENE", EventType.DELETION));
    }

    @NotNull
    private static CopyNumberExtractor createTestExtractor(@NotNull DriverInconsistencyMode mode) {
        return new CopyNumberExtractor(new GeneChecker(Sets.newHashSet("PTEN", "AKT1", "KRAS")),
                DriverGeneTestFactory.createDriverGenes("KRAS", "AKT1"),
                mode);
    }
}