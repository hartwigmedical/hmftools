package com.hartwig.hmftools.common.drivercatalog;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelFactoryTest;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CNADriversTest {

    private final DriverGenePanel genePanel = DriverGenePanelFactoryTest.testGenePanel();

    @Test
    public void testDeletionsInGermlineAsStillReportable() {
        GeneCopyNumber del = createTestCopyNumberBuilder("APC").germlineHet2HomRegions(1).germlineHomRegions(1).build();
        List<DriverCatalog> drivers = new CNADrivers(genePanel).deletions(Lists.newArrayList(del));
        assertEquals(1, drivers.size());
        assertEquals("APC", drivers.get(0).gene());
    }

    @NotNull
    private static ImmutableGeneCopyNumber.Builder createTestCopyNumberBuilder(@NotNull String gene) {
        return ImmutableGeneCopyNumber.builder()
                .start(1)
                .end(2)
                .gene(gene)
                .chromosome("1")
                .chromosomeBand("band")
                .minRegionStart(0)
                .minRegionStartSupport(SegmentSupport.NONE)
                .minRegionEnd(0)
                .minRegionEndSupport(SegmentSupport.NONE)
                .minRegionMethod(CopyNumberMethod.UNKNOWN)
                .minRegions(1)
                .germlineHet2HomRegions(0)
                .germlineHomRegions(0)
                .somaticRegions(1)
                .minCopyNumber(0.1)
                .maxCopyNumber(0.1)
                .transcriptID("trans")
                .transcriptVersion(0)
                .minMinorAlleleCopyNumber(0);
    }
}
