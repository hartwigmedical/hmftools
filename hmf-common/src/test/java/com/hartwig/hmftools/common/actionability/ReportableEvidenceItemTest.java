package com.hartwig.hmftools.common.actionability;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Ignore;
import org.junit.Test;

public class ReportableEvidenceItemTest {

    @Test
    public void canFilterEvidenceItems() {
        List<EvidenceItem> evidenceItems = Lists.newArrayList();

        evidenceItems.add(testEvidenceBuilder().event("TP53 p.Val600Glu").drug("Tamoxifen").build());
        evidenceItems.add(testEvidenceBuilder().event("BRAF p.Val600Glu").drug("Cobimetinib + Vemurafenib").build());

        assertEquals(1, ReportableEvidenceItem.extractAllReportableEvidenceItems(evidenceItems).size());
    }

    @Test
    @Ignore
    public void canFilterEvidenceItems2() {
        List<EvidenceItem> evidenceItems = Lists.newArrayList();

        evidenceItems.add(testEvidenceBuilder().event("TP53 Deletion").drug("Nivolumab").build());
        evidenceItems.add(testEvidenceBuilder().event("BRAF p.Val600Glu").drug("Tamoxifen").build());

        // TODO This should result 2 evidence items!
        assertEquals(2, ReportableEvidenceItem.extractAllReportableEvidenceItems(evidenceItems).size());
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder testEvidenceBuilder() {
        return ImmutableEvidenceItem.builder()
                .isOnLabel(true)
                .drugsType(Strings.EMPTY)
                .cancerType(Strings.EMPTY)
                .response(Strings.EMPTY)
                .level(EvidenceLevel.LEVEL_A)
                .reference(Strings.EMPTY)
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.SPECIFIC);
    }

}