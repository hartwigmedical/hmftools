package com.hartwig.hmftools.common.drivercatalog;

import static com.hartwig.hmftools.common.purple.Gender.FEMALE;
import static com.hartwig.hmftools.common.purple.Gender.MALE;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelFactoryTest;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberTestFactory;
import com.hartwig.hmftools.common.purple.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CopyNumberDriversTest
{
    private final DriverGenePanel genePanel = DriverGenePanelFactoryTest.testGenePanel();

    @Test
    public void testOncoAmplificationWithoutDnds()
    {
        DriverGene onco = genePanel.amplificationTargets()
                .stream()
                .filter(x -> !x.reportSomatic() && x.likelihoodType() == DriverCategory.ONCO)
                .findFirst()
                .get();

        GeneCopyNumber oncoAmp = createTestCopyNumberBuilder(onco.gene()).minCopyNumber(7).maxCopyNumber(7).build();

        List<DriverCatalog> drivers = AmplificationDrivers.findAmplifications(
                Sets.newHashSet(PurpleQCStatus.PASS), FEMALE, genePanel, 2, Lists.newArrayList(oncoAmp), false);

        assertEquals(oncoAmp.geneName(), drivers.get(0).gene());
        assertEquals(DriverCategory.ONCO, drivers.get(0).category());

        // test amplification on chrX for MALE
        GeneCopyNumber chrXAmp = createTestCopyNumberBuilder("AR")
                .chromosome(HumanChromosome._X.toString()).minCopyNumber(7).maxCopyNumber(7).build();

        drivers = AmplificationDrivers.findAmplifications(
                Sets.newHashSet(PurpleQCStatus.PASS), MALE, genePanel, 2, Lists.newArrayList(chrXAmp), false);

        assertEquals(chrXAmp.geneName(), drivers.get(0).gene());
        assertEquals(DriverCategory.ONCO, drivers.get(0).category());
    }

    @Test
    public void testPartialAmp()
    {
        List<DriverGene> driverGenes = Lists.newArrayList(genePanel.amplificationTargets());
        GeneCopyNumber partialAmp = createTestCopyNumberBuilder(driverGenes.get(0).gene()).minCopyNumber(0.1).maxCopyNumber(7).build();
        GeneCopyNumber fullAmp = createTestCopyNumberBuilder(driverGenes.get(1).gene()).minCopyNumber(7).maxCopyNumber(7).build();

        List<DriverCatalog> drivers = AmplificationDrivers.findAmplifications(
                Sets.newHashSet(PurpleQCStatus.PASS), FEMALE, genePanel,
                2, Lists.newArrayList(partialAmp, fullAmp), false);
        assertEquals(2, drivers.size());

        assertEquals(partialAmp.geneName(), drivers.get(0).gene());
        assertEquals(DriverType.PARTIAL_AMP, drivers.get(0).driver());

        assertEquals(fullAmp.geneName(), drivers.get(1).gene());
        assertEquals(DriverType.AMP, drivers.get(1).driver());
    }

    @Test
    public void testAmpWithWarn()
    {
        String gene = Lists.newArrayList(genePanel.amplificationTargets()).get(0).gene();
        GeneCopyNumber ampNoSupport = createTestCopyNumberBuilder(gene).geneName(gene)
                .minCopyNumber(100)
                .minRegionStartSupport(SegmentSupport.NONE)
                .minRegionEndSupport(SegmentSupport.NONE)
                .build();

        GeneCopyNumber ampStartSupport =
                ImmutableGeneCopyNumber.builder().from(ampNoSupport).minRegionStartSupport(SegmentSupport.BND).build();
        GeneCopyNumber ampEndSupport = ImmutableGeneCopyNumber.builder().from(ampNoSupport).minRegionEndSupport(SegmentSupport.BND).build();
        GeneCopyNumber ampBothSupport = ImmutableGeneCopyNumber.builder()
                .from(ampNoSupport)
                .minRegionStartSupport(SegmentSupport.BND)
                .minRegionEndSupport(SegmentSupport.BND)
                .build();

        Set<PurpleQCStatus> warnDeletedGenes = Sets.newHashSet(PurpleQCStatus.WARN_DELETED_GENES);
        assertEquals(1, AmplificationDrivers.findAmplifications(
                warnDeletedGenes, FEMALE, genePanel,1, Collections.singletonList(ampNoSupport), false).size());

        Set<PurpleQCStatus> warnCopyNumber = Sets.newHashSet(PurpleQCStatus.WARN_HIGH_COPY_NUMBER_NOISE);
        assertEquals(0, AmplificationDrivers.findAmplifications(
                warnCopyNumber, FEMALE, genePanel,1, Collections.singletonList(ampNoSupport), false).size());

        assertEquals(1, AmplificationDrivers.findAmplifications(
                warnCopyNumber, FEMALE, genePanel,1, Collections.singletonList(ampStartSupport), false).size());

        assertEquals(1, AmplificationDrivers.findAmplifications(
                warnCopyNumber, FEMALE, genePanel,1, Collections.singletonList(ampEndSupport), false).size());

        assertEquals(1, AmplificationDrivers.findAmplifications(
                warnCopyNumber, FEMALE, genePanel,1, Collections.singletonList(ampBothSupport), false).size());
    }

    @Test
    public void testDelWithWarn()
    {
        String gene = Lists.newArrayList(genePanel.deletionTargets()).get(0).gene();
        GeneCopyNumber longDelNoSupport = createTestCopyNumberBuilder(gene).geneName(gene)
                .minRegionStart(1)
                .minRegionEnd(100_000_000)
                .minCopyNumber(0.01)
                .minRegionStartSupport(SegmentSupport.NONE)
                .minRegionEndSupport(SegmentSupport.NONE)
                .build();

        GeneCopyNumber longDelStartSupport = ImmutableGeneCopyNumber.builder()
                .from(longDelNoSupport)
                .minRegionStartSupport(SegmentSupport.BND)
                .minRegionEndSupport(SegmentSupport.TELOMERE)
                .build();
        GeneCopyNumber longDelEndSupport = ImmutableGeneCopyNumber.builder()
                .from(longDelNoSupport)
                .minRegionStartSupport(SegmentSupport.CENTROMERE)
                .minRegionStartSupport(SegmentSupport.BND)
                .build();

        GeneCopyNumber shortDelStartSupport = ImmutableGeneCopyNumber.builder()
                .from(longDelNoSupport)
                .minRegionStartSupport(SegmentSupport.BND)
                .minRegionEndSupport(SegmentSupport.TELOMERE)
                .minRegionStart(1)
                .minRegionEnd(100)
                .build();

        GeneCopyNumber shortDelEndSupport = ImmutableGeneCopyNumber.builder()
                .from(longDelNoSupport)
                .minRegionStartSupport(SegmentSupport.CENTROMERE)
                .minRegionEndSupport(SegmentSupport.BND)
                .minRegionStart(1)
                .minRegionEnd(100)
                .build();

        GeneCopyNumber longDelBothSupport = ImmutableGeneCopyNumber.builder()
                .from(longDelNoSupport)
                .minRegionStartSupport(SegmentSupport.BND)
                .minRegionEndSupport(SegmentSupport.BND)
                .build();

        Set<PurpleQCStatus> noWarn = Sets.newHashSet();
        assertEquals(1, new DeletionDrivers(noWarn, genePanel).deletions(Collections.singletonList(longDelNoSupport), false).size());

        Set<PurpleQCStatus> warnDeletedGenes = Sets.newHashSet(PurpleQCStatus.WARN_DELETED_GENES);

        assertEquals(0, new DeletionDrivers(warnDeletedGenes, genePanel).deletions(
                Collections.singletonList(longDelNoSupport), false).size());

        assertEquals(0, new DeletionDrivers(warnDeletedGenes, genePanel).deletions(
                Collections.singletonList(longDelStartSupport), false).size());

        assertEquals(0, new DeletionDrivers(warnDeletedGenes, genePanel).deletions(
                Collections.singletonList(longDelEndSupport), false).size());

        assertEquals(1, new DeletionDrivers(warnDeletedGenes, genePanel).deletions(
                Collections.singletonList(shortDelStartSupport), false).size());

        assertEquals(1, new DeletionDrivers(warnDeletedGenes, genePanel).deletions(
                Collections.singletonList(shortDelEndSupport), false).size());

        assertEquals(1, new DeletionDrivers(warnDeletedGenes, genePanel).deletions(
                Collections.singletonList(longDelBothSupport), false).size());

        Set<PurpleQCStatus> warnCopyNumber = Sets.newHashSet(PurpleQCStatus.WARN_HIGH_COPY_NUMBER_NOISE);
        assertEquals(0, new DeletionDrivers(warnCopyNumber, genePanel).deletions(
                Collections.singletonList(longDelNoSupport), false).size());

        assertEquals(0, new DeletionDrivers(warnCopyNumber, genePanel).deletions(
                Collections.singletonList(longDelStartSupport), false).size());

        assertEquals(0, new DeletionDrivers(warnCopyNumber, genePanel).deletions(
                Collections.singletonList(longDelEndSupport), false).size());

        assertEquals(1, new DeletionDrivers(warnCopyNumber, genePanel).deletions(
                Collections.singletonList(shortDelStartSupport), false).size());

        assertEquals(1, new DeletionDrivers(warnCopyNumber, genePanel).deletions(
                Collections.singletonList(shortDelEndSupport), false).size());

        assertEquals(1, new DeletionDrivers(warnCopyNumber, genePanel).deletions(
                Collections.singletonList(longDelBothSupport), false).size());
    }

    @NotNull
    private static ImmutableGeneCopyNumber.Builder createTestCopyNumberBuilder(@NotNull String gene)
    {
        return GeneCopyNumberTestFactory.builder().geneName(gene);
    }
}
