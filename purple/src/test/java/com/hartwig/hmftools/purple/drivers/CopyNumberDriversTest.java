package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.common.purple.Gender.FEMALE;
import static com.hartwig.hmftools.common.purple.Gender.MALE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGenePanel;
import com.hartwig.hmftools.common.driver.panel.DriverGenePanelFactoryTest;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberTestFactory;
import com.hartwig.hmftools.common.purple.SegmentSupport;

import org.apache.logging.log4j.util.Strings;
import org.junit.Assert;
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

        GeneCopyNumber oncoAmp = GeneCopyNumberTestFactory.createGeneCopyNumber(onco.gene(), 7, 7);

        List<DriverCatalog> drivers = AmplificationDrivers.findAmplifications(
                Sets.newHashSet(PurpleQCStatus.PASS), FEMALE, genePanel, 2, Lists.newArrayList(oncoAmp), false);

        assertEquals(oncoAmp.geneName(), drivers.get(0).gene());
        assertEquals(DriverCategory.ONCO, drivers.get(0).category());

        // test amplification on chrX for MALE
        GeneCopyNumber chrXAmp = GeneCopyNumberTestFactory.createGeneCopyNumber(
                HumanChromosome._X.toString(), "AR", 7, 7);

        drivers = AmplificationDrivers.findAmplifications(
                Sets.newHashSet(PurpleQCStatus.PASS), MALE, genePanel, 2, Lists.newArrayList(chrXAmp), false);

        assertEquals(chrXAmp.geneName(), drivers.get(0).gene());
        assertEquals(DriverCategory.ONCO, drivers.get(0).category());
    }

    @Test
    public void testPartialAmp()
    {
        List<DriverGene> driverGenes = Lists.newArrayList(genePanel.amplificationTargets());

        GeneCopyNumber partialAmp = GeneCopyNumberTestFactory.createGeneCopyNumber(
                driverGenes.get(0).gene(), 0.1, 7);

        GeneCopyNumber fullAmp = GeneCopyNumberTestFactory.createGeneCopyNumber(
                driverGenes.get(1).gene(), 7, 7);

        List<DriverCatalog> drivers = AmplificationDrivers.findAmplifications(
                Sets.newHashSet(PurpleQCStatus.PASS), FEMALE, genePanel,
                2, Lists.newArrayList(partialAmp, fullAmp), false);
        assertEquals(2, drivers.size());

        assertEquals(partialAmp.geneName(), drivers.get(0).gene());
        Assert.assertEquals(DriverType.PARTIAL_AMP, drivers.get(0).driver());

        assertEquals(fullAmp.geneName(), drivers.get(1).gene());
        assertEquals(DriverType.AMP, drivers.get(1).driver());
    }

    @Test
    public void testAmpWithWarn()
    {
        String gene = Lists.newArrayList(genePanel.amplificationTargets()).get(0).gene();

        GeneCopyNumber ampNoSupport = createGeneCopyNumber(
                gene, 100, 1, 10000, SegmentSupport.NONE, SegmentSupport.NONE);

        GeneCopyNumber ampStartSupport = createGeneCopyNumber(
                gene, 100, 1, 10000, SegmentSupport.BND, SegmentSupport.NONE);

        GeneCopyNumber ampEndSupport = createGeneCopyNumber(
                gene, 100, 1, 10000, SegmentSupport.NONE, SegmentSupport.BND);

        GeneCopyNumber ampBothSupport = createGeneCopyNumber(
                gene, 100, 1, 10000, SegmentSupport.BND, SegmentSupport.BND);

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

        GeneCopyNumber longDelNoSupport = createGeneCopyNumber(
                gene, 0.01, 1, 100_000_000, SegmentSupport.NONE, SegmentSupport.NONE);

        GeneCopyNumber longDelStartSupport = createGeneCopyNumber(
                gene, 0.01, 1, 100_000_000, SegmentSupport.BND, SegmentSupport.TELOMERE);

        GeneCopyNumber longDelEndSupport = createGeneCopyNumber(
                gene, 0.01, 1, 100_000_000, SegmentSupport.CENTROMERE, SegmentSupport.BND);

        GeneCopyNumber longDelBothSupport = createGeneCopyNumber(
                gene, 0.01, 1, 100_000_000, SegmentSupport.BND, SegmentSupport.BND);

        GeneCopyNumber shortDelNoSupport = createGeneCopyNumber(
                gene, 0.01, 1, 100, SegmentSupport.NONE, SegmentSupport.NONE);

        GeneCopyNumber shortDelStartSupport = createGeneCopyNumber(
                gene, 0.01, 1, 100, SegmentSupport.BND, SegmentSupport.TELOMERE);

        GeneCopyNumber shortDelEndSupport = createGeneCopyNumber(
                gene, 0.01, 1, 100, SegmentSupport.CENTROMERE, SegmentSupport.BND);

        GeneCopyNumber shortDelBothSupport = createGeneCopyNumber(
                gene, 0.01, 1, 100, SegmentSupport.BND, SegmentSupport.BND);

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
                Collections.singletonList(longDelBothSupport), false).size());

        assertEquals(0, new DeletionDrivers(warnDeletedGenes, genePanel).deletions(
                Collections.singletonList(longDelNoSupport), false).size());

        assertEquals(0, new DeletionDrivers(warnDeletedGenes, genePanel).deletions(
                Collections.singletonList(shortDelNoSupport), false).size());

        assertEquals(1, new DeletionDrivers(warnDeletedGenes, genePanel).deletions(
                Collections.singletonList(shortDelEndSupport), false).size());

        assertEquals(1, new DeletionDrivers(warnDeletedGenes, genePanel).deletions(
                Collections.singletonList(shortDelBothSupport), false).size());

        // Check effect of panel mode.
        assertEquals(0, new DeletionDrivers(warnDeletedGenes, genePanel).deletions(
                Collections.singletonList(longDelNoSupport), true).size());

        assertEquals(0, new DeletionDrivers(warnDeletedGenes, genePanel).deletions(
                Collections.singletonList(longDelStartSupport), true).size());

        assertEquals(0, new DeletionDrivers(warnDeletedGenes, genePanel).deletions(
                Collections.singletonList(longDelEndSupport), true).size());

        assertEquals(1, new DeletionDrivers(warnDeletedGenes, genePanel).deletions(
                Collections.singletonList(longDelBothSupport), true).size());

        assertEquals(0, new DeletionDrivers(warnDeletedGenes, genePanel).deletions(
                Collections.singletonList(shortDelNoSupport), true).size());

        assertEquals(0, new DeletionDrivers(warnDeletedGenes, genePanel).deletions(
                Collections.singletonList(shortDelStartSupport), true).size());

        assertEquals(0, new DeletionDrivers(warnDeletedGenes, genePanel).deletions(
                Collections.singletonList(shortDelEndSupport), true).size());

        assertEquals(1, new DeletionDrivers(warnDeletedGenes, genePanel).deletions(
                Collections.singletonList(shortDelBothSupport), true).size());

        Set<PurpleQCStatus> warnCopyNumber = Sets.newHashSet(PurpleQCStatus.WARN_HIGH_COPY_NUMBER_NOISE);
        assertEquals(0, new DeletionDrivers(warnCopyNumber, genePanel).deletions(
                Collections.singletonList(longDelNoSupport), false).size());

        assertEquals(0, new DeletionDrivers(warnCopyNumber, genePanel).deletions(
                Collections.singletonList(longDelStartSupport), false).size());

        assertEquals(0, new DeletionDrivers(warnCopyNumber, genePanel).deletions(
                Collections.singletonList(longDelEndSupport), false).size());

        assertEquals(1, new DeletionDrivers(warnCopyNumber, genePanel).deletions(
                Collections.singletonList(longDelBothSupport), false).size());

        assertEquals(0, new DeletionDrivers(warnDeletedGenes, genePanel).deletions(
                Collections.singletonList(shortDelNoSupport), false).size());

        assertEquals(1, new DeletionDrivers(warnCopyNumber, genePanel).deletions(
                Collections.singletonList(shortDelStartSupport), false).size());

        assertEquals(1, new DeletionDrivers(warnCopyNumber, genePanel).deletions(
                Collections.singletonList(shortDelEndSupport), false).size());

        assertEquals(1, new DeletionDrivers(warnCopyNumber, genePanel).deletions(
                Collections.singletonList(shortDelBothSupport), false).size());
    }

    private static GeneCopyNumber createGeneCopyNumber(
            final String gene, double minCopyNumber, int minRegionStart, int minRegionEnd,
            final SegmentSupport minRegionStartSupport, final SegmentSupport minRegionEndSupport)
    {
        return new GeneCopyNumber(CHR_1, 0, 0, gene, Strings.EMPTY, true,
                Strings.EMPTY, 0 ,minCopyNumber, 0, 1, 1,
                minRegionStart, minRegionEnd, 0, 1.0,
                minRegionStartSupport, minRegionEndSupport, CopyNumberMethod.UNKNOWN);
    }
}
