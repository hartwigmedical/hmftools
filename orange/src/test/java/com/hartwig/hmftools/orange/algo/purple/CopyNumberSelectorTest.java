package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneTestFactory;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberTestFactory;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;

import org.junit.Test;

public class CopyNumberSelectorTest
{
    @Test
    public void canSelectNearReportableGains()
    {
        DriverGene driver1 = DriverGeneTestFactory.builder().gene("driver 1").reportAmplification(true).build();
        DriverGene driver2 = DriverGeneTestFactory.builder().gene("driver 2").reportAmplification(true).build();
        DriverGene driver3 = DriverGeneTestFactory.builder().gene("driver 3").reportAmplification(false).build();
        List<DriverGene> driverGenes = Lists.newArrayList(driver1, driver2, driver3);

        GeneCopyNumber match = GeneCopyNumberTestFactory.createGeneCopyNumber(driver1.gene(), 11D, 11D);
        GeneCopyNumber tooLowCN = GeneCopyNumberTestFactory.createGeneCopyNumber(driver2.gene(), 5D, 5D);
        GeneCopyNumber tooHighCN =
                GeneCopyNumberTestFactory.createGeneCopyNumber(driver2.gene(), 15D, 15D);
        GeneCopyNumber noReportAmp =
                GeneCopyNumberTestFactory.createGeneCopyNumber(driver3.gene(), 11D, 11D);
        List<GeneCopyNumber> allGeneCopyNumbers = Lists.newArrayList(match, tooLowCN, tooHighCN, noReportAmp);

        List<PurpleGainDeletion> nearReportableGains =
                CopyNumberSelector.selectNearReportableSomaticGains(allGeneCopyNumbers, 4D, Lists.newArrayList(), driverGenes);

        assertEquals(1, nearReportableGains.size());
        assertEquals(driver1.gene(), nearReportableGains.get(0).gene());
    }

    @Test
    public void warnIfNearReportableIsReported()
    {
        DriverGene driver = DriverGeneTestFactory.builder().gene("driver 1").reportAmplification(true).build();

        GeneCopyNumber match = GeneCopyNumberTestFactory.createGeneCopyNumber(driver.gene(), 11D, 0);

        PurpleGainDeletion reportable = TestPurpleGainDeletionFactory.builder().gene(match.geneName()).build();

        assertNotNull(CopyNumberSelector.selectNearReportableSomaticGains(Lists.newArrayList(match),
                4D,
                Lists.newArrayList(reportable),
                Lists.newArrayList(driver)));
    }
}