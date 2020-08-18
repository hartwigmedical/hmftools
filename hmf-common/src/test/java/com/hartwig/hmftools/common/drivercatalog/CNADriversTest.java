package com.hartwig.hmftools.common.drivercatalog;

import static org.junit.Assert.assertEquals;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelFactory;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CNADriversTest {

    private final DriverGenePanel genePanel = new DriverGenePanelFactory().create();

    @Test
    public void testDeletionsInGermlineAsStillReportable() {
        GeneCopyNumber del = createTestCopyNumberBuilder("APC").germlineHet2HomRegions(1).germlineHomRegions(1).build();
        List<DriverCatalog> drivers = new CNADrivers(genePanel).deletions(Lists.newArrayList(del));
        assertEquals(1, drivers.size());
        assertEquals("APC", drivers.get(0).gene());
    }

    @Test
    public void testChromosomeBand() {
        GeneCopyNumber mapped = createTestCopyNumberBuilder("AC093642.5").build();
        GeneCopyNumber unmapped = createTestCopyNumberBuilder("APC").build();

        List<DriverCatalog> drivers = new CNADrivers(genePanel).deletions(Lists.newArrayList(mapped, unmapped));
        assertEquals(2, drivers.size());
        assertEquals("q telomere", drivers.get(0).chromosomeBand());
        assertEquals("band", drivers.get(1).chromosomeBand());
    }

    @Test
    public void testAmpTargets() {
        Set<String> oldTargets = amplificationTargets();
        Set<String> newTargets = genePanel.amplificationTargets();
        assertEquals(oldTargets.size(), newTargets.size());
        assertEquals(0, Sets.difference(oldTargets, newTargets).size());
    }

    @Test
    public void testDelTargets() {
        Set<String> oldTargets = deletionTargets().keySet();
        Set<String> newTargets = genePanel.deletionTargets();
        assertEquals(oldTargets.size(), newTargets.size());
        assertEquals(0, Sets.difference(oldTargets, newTargets).size());
    }

    @Test
    public void testDelBandTargets() {
        Map<String, String> oldTargets = deletionTargets();
        Map<String, String> newTargets = genePanel.deletionBandMap();
        for (Map.Entry<String, String> entry : newTargets.entrySet()) {
            assertEquals(oldTargets.get(entry.getKey()), entry.getValue());
        }
    }

    @NotNull
    private static Set<String> amplificationTargets() {
        final InputStream inputStream = DndsDriverGeneLikelihoodSupplier.class.getResourceAsStream("/cna/AmplificationTargets.tsv");
        return new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.toSet());
    }

    @NotNull
    private static Map<String, String> deletionTargets() {
        final Map<String, String> result = Maps.newHashMap();
        final InputStream inputStream = DndsDriverGeneLikelihoodSupplier.class.getResourceAsStream("/cna/DeletionTargets.tsv");
        new BufferedReader(new InputStreamReader(inputStream)).lines().forEach(line -> {
            final String[] values = line.split("\t");
            result.put(values[0], values[1]);
        });

        return result;
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
