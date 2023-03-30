package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneTestFactory;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberTestFactory;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGainLoss;
import com.hartwig.hmftools.datamodel.purple.PurpleGainLoss;

import org.junit.Test;

public class CopyNumberSelectorTest {

    @Test
    public void canSelectNearReportableGains() {
        DriverGene driver1 = DriverGeneTestFactory.builder().gene("driver 1").reportAmplification(true).build();
        DriverGene driver2 = DriverGeneTestFactory.builder().gene("driver 2").reportAmplification(true).build();
        DriverGene driver3 = DriverGeneTestFactory.builder().gene("driver 3").reportAmplification(false).build();
        List<DriverGene> driverGenes = Lists.newArrayList(driver1, driver2, driver3);

        GeneCopyNumber match = GeneCopyNumberTestFactory.builder().geneName(driver1.gene()).minCopyNumber(11D).maxCopyNumber(11D).build();
        GeneCopyNumber tooLowCN = GeneCopyNumberTestFactory.builder().geneName(driver2.gene()).minCopyNumber(5D).maxCopyNumber(5D).build();
        GeneCopyNumber tooHighCN =
                GeneCopyNumberTestFactory.builder().geneName(driver2.gene()).minCopyNumber(15D).maxCopyNumber(15D).build();
        GeneCopyNumber noReportAmp =
                GeneCopyNumberTestFactory.builder().geneName(driver3.gene()).minCopyNumber(11D).maxCopyNumber(11D).build();
        List<GeneCopyNumber> allGeneCopyNumbers = Lists.newArrayList(match, tooLowCN, tooHighCN, noReportAmp);

        List<PurpleGainLoss> nearReportableGains =
                CopyNumberSelector.selectNearReportableSomaticGains(allGeneCopyNumbers, 4D, Lists.newArrayList(), driverGenes);

        assertEquals(1, nearReportableGains.size());
        assertEquals(driver1.gene(), nearReportableGains.get(0).gene());
    }

    @Test
    public void warnIfNearReportableIsReported() {
        DriverGene driver = DriverGeneTestFactory.builder().gene("driver 1").reportAmplification(true).build();

        GeneCopyNumber match = GeneCopyNumberTestFactory.builder().geneName(driver.gene()).minCopyNumber(11D).build();

        PurpleGainLoss reportable = TestPurpleGainLossFactory.builder().gene(match.geneName()).build();

        assertNotNull(CopyNumberSelector.selectNearReportableSomaticGains(Lists.newArrayList(match),
                4D,
                Lists.newArrayList(reportable),
                Lists.newArrayList(driver)));
    }

    @Test
    public void canSelectPotentiallyInterestingGains() {
        PurpleGainLoss gain = TestPurpleGainLossFactory.builder()
                .chromosome("chr1")
                .chromosomeBand("band 1")
                .gene("gene 1")
                .interpretation(CopyNumberInterpretation.FULL_GAIN)
                .minCopies(12.1)
                .maxCopies(12.1)
                .build();
        PurpleGainLoss lowerGainSameBand = TestPurpleGainLossFactory.builder()
                .chromosome("chr1")
                .chromosomeBand("band 1")
                .gene("gene 2")
                .interpretation(CopyNumberInterpretation.FULL_GAIN)
                .minCopies(11.1)
                .maxCopies(11.1)
                .build();
        PurpleGainLoss lowestGainOtherBand = TestPurpleGainLossFactory.builder()
                .chromosome("chr1")
                .chromosomeBand("band 2")
                .gene("gene 3")
                .interpretation(CopyNumberInterpretation.FULL_GAIN)
                .minCopies(10.1)
                .maxCopies(10.1)
                .build();
        PurpleGainLoss partialGain = TestPurpleGainLossFactory.builder()
                .chromosome("chr1")
                .chromosomeBand("band 3")
                .gene("gene 4")
                .interpretation(CopyNumberInterpretation.PARTIAL_GAIN)
                .minCopies(4.1)
                .maxCopies(15.1)
                .build();

        List<PurpleGainLoss> allGains = Lists.newArrayList(gain, lowerGainSameBand, lowestGainOtherBand, partialGain);

        List<PurpleGainLoss> interesting = CopyNumberSelector.selectInterestingUnreportedGainsLosses(allGains, Lists.newArrayList());

        assertEquals(2, interesting.size());
        assertTrue(interesting.contains(gain));
        assertTrue(interesting.contains(lowestGainOtherBand));
    }

    @Test
    public void canSelectPotentiallyInterestingLosses() {
        ImmutablePurpleGainLoss.Builder lossBuilder =
                TestPurpleGainLossFactory.builder().interpretation(CopyNumberInterpretation.FULL_LOSS).minCopies(0).maxCopies(0);

        PurpleGainLoss interestingLoss = lossBuilder.chromosome("chr2").chromosomeBand("band 1").gene("gene 2").build();
        PurpleGainLoss lossWithReportable = lossBuilder.chromosome("chr2").chromosomeBand("band 2").gene("gene 3").build();
        PurpleGainLoss reportableLossSameBand = lossBuilder.chromosome("chr2").chromosomeBand("band 2").gene("gene 4").build();
        PurpleGainLoss lossAllosome = lossBuilder.chromosome("chrX").chromosomeBand("band 1").gene("gene 5").build();
        PurpleGainLoss reportableLossOtherBand = lossBuilder.chromosome("chr2").chromosomeBand("band 3").gene("gene 6").build();

        List<PurpleGainLoss> allLosses =
                Lists.newArrayList(interestingLoss, lossWithReportable, reportableLossOtherBand, lossAllosome, reportableLossOtherBand);

        List<PurpleGainLoss> interesting = CopyNumberSelector.selectInterestingUnreportedGainsLosses(allLosses,
                Lists.newArrayList(reportableLossOtherBand, reportableLossSameBand));

        assertEquals(1, interesting.size());
        assertTrue(interesting.contains(interestingLoss));
    }
}