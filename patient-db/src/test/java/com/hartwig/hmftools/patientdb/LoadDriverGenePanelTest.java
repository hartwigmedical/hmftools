package com.hartwig.hmftools.patientdb;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.driver.panel.ImmutableDriverGene;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.DrivergenepanelRecord;

import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class LoadDriverGenePanelTest extends DatabaseTestBase
{
    @Test
    public void canLoadDriverGenePanel()
    {
        List<DriverGene> driverGenes = List.of(
                createTestGene("BRCA1"),
                createTestGene("BRCA2")
        );

        databaseAccess.writeGenePanel(driverGenes);

        List<DrivergenepanelRecord> driverGeneRecords = fetchTable(Tables.DRIVERGENEPANEL, DrivergenepanelRecord.class);

        assertEquals(2, driverGeneRecords.size());
    }

    private static DriverGene createTestGene(String geneName)
    {
        return ImmutableDriverGene.builder()
                .gene(geneName)
                .reportMissenseAndInframe(false)
                .reportNonsenseAndFrameshift(false)
                .reportSplice(false)
                .reportDeletion(false)
                .reportHetDeletion(false)
                .hetDeletionThreshold(0)
                .reportDisruption(false)
                .reportAmplification(false)
                .amplificationRatio(0)
                .reportSomaticHotspot(false)
                .reportGermlineVariant(DriverGeneGermlineReporting.NONE)
                .reportGermlineHotspot(DriverGeneGermlineReporting.NONE)
                .likelihoodType(DriverCategory.TSG)
                .reportGermlineDeletion(DriverGeneGermlineReporting.NONE)
                .reportGermlineDisruption(DriverGeneGermlineReporting.NONE)
                .additionalReportedTranscripts(List.of())
                .reportPGX(false)
                .build();
    }
}
