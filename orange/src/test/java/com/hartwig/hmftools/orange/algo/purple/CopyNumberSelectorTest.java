package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneTestFactory;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberTestFactory;
import com.hartwig.hmftools.common.purple.interpretation.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.interpretation.GainLoss;
import com.hartwig.hmftools.common.purple.interpretation.GainLossTestFactory;

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

        List<GainLoss> nearReportableGains =
                CopyNumberSelector.selectNearReportableSomaticGains(allGeneCopyNumbers, 4D, Lists.newArrayList(), driverGenes);

        assertEquals(1, nearReportableGains.size());
        assertEquals(driver1.gene(), nearReportableGains.get(0).gene());
    }

    @Test
    public void warnIfNearReportableIsReported() {
        DriverGene driver = DriverGeneTestFactory.builder().gene("driver 1").reportAmplification(true).build();

        GeneCopyNumber match = GeneCopyNumberTestFactory.builder().geneName(driver.gene()).minCopyNumber(11D).build();

        GainLoss reportable = GainLossTestFactory.builder().gene(match.geneName()).build();

        assertNotNull(CopyNumberSelector.selectNearReportableSomaticGains(Lists.newArrayList(match),
                4D,
                Lists.newArrayList(reportable),
                Lists.newArrayList(driver)));
    }

    @Test
    public void canSelectPotentiallyInterestingGainsLosses() {
        GainLoss gain = GainLossTestFactory.builder()
                .chromosome("chr1")
                .chromosomeBand("band 1")
                .gene("gene 1")
                .interpretation(CopyNumberInterpretation.FULL_GAIN)
                .minCopies(12)
                .maxCopies(12)
                .build();
        GainLoss lowerGainSameBand = GainLossTestFactory.builder()
                .chromosome("chr1")
                .chromosomeBand("band 1")
                .gene("gene 2")
                .interpretation(CopyNumberInterpretation.FULL_GAIN)
                .minCopies(11)
                .maxCopies(11)
                .build();
        GainLoss lowestGainOtherBand = GainLossTestFactory.builder()
                .chromosome("chr1")
                .chromosomeBand("band 2")
                .gene("gene 3")
                .interpretation(CopyNumberInterpretation.FULL_GAIN)
                .minCopies(10)
                .maxCopies(10)
                .build();
        GainLoss partialGain = GainLossTestFactory.builder()
                .chromosome("chr1")
                .chromosomeBand("band 3")
                .gene("gene 4")
                .interpretation(CopyNumberInterpretation.PARTIAL_GAIN)
                .minCopies(4)
                .maxCopies(15)
                .build();

        GainLoss loss = GainLossTestFactory.builder()
                .chromosome("chr2")
                .chromosomeBand("band 1")
                .gene("gene 5")
                .interpretation(CopyNumberInterpretation.FULL_LOSS)
                .minCopies(0)
                .maxCopies(0)
                .build();
        GainLoss otherLoss = GainLossTestFactory.builder()
                .chromosome("chr2")
                .chromosomeBand("band 2")
                .gene("gene 6")
                .interpretation(CopyNumberInterpretation.FULL_LOSS)
                .minCopies(0)
                .maxCopies(0)
                .build();
        GainLoss reportableLossSameBand = GainLossTestFactory.builder()
                .chromosome("chr2")
                .chromosomeBand("band 2")
                .gene("gene 7")
                .interpretation(CopyNumberInterpretation.FULL_LOSS)
                .minCopies(0)
                .maxCopies(0)
                .build();
        GainLoss lossAllosome = GainLossTestFactory.builder()
                .chromosome("chrX")
                .chromosomeBand("band 1")
                .gene("gene 8")
                .interpretation(CopyNumberInterpretation.FULL_LOSS)
                .minCopies(0)
                .maxCopies(0)
                .build();
        GainLoss reportableLossOtherBand = GainLossTestFactory.builder()
                .chromosome("chr2")
                .chromosomeBand("band 3")
                .gene("gene 9")
                .interpretation(CopyNumberInterpretation.FULL_LOSS)
                .minCopies(0)
                .maxCopies(0)
                .build();

        List<GainLoss> allGainsLosses = Lists.newArrayList(gain,
                lowerGainSameBand,
                lowestGainOtherBand,
                partialGain,
                loss,
                otherLoss,
                reportableLossOtherBand,
                lossAllosome,
                reportableLossOtherBand);

        List<GainLoss> interesting = CopyNumberSelector.selectInterestingUnreportedGainsLosses(allGainsLosses,
                Lists.newArrayList(reportableLossOtherBand, reportableLossSameBand));

        assertEquals(3, interesting.size());
        assertTrue(interesting.contains(gain));
        assertTrue(interesting.contains(lowestGainOtherBand));
        assertTrue(interesting.contains(loss));
    }
}