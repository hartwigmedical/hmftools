package com.hartwig.hmftools.patientdb;


import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.virus.ImmutableVirusBreakend;
import com.hartwig.hmftools.common.virus.VirusBreakend;
import com.hartwig.hmftools.common.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.VirusbreakendRecord;

import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class LoadVirusBreakendDataTest extends DatabaseTestBase
{
    @Test
    public void canWriteVirusBreakendData()
    {
        VirusBreakend virusBreakend = ImmutableVirusBreakend.builder()
                .taxidGenus(333750)
                .nameGenus("Alphapapillomavirus")
                .readsGenusTree(14884)
                .taxidSpecies(337041)
                .nameSpecies("Alphapapillomavirus 9")
                .readsSpeciesTree(14883)
                .taxidAssigned(333760)
                .nameAssigned("Human papillomavirus 16")
                .readsAssignedTree(14882)
                .readsAssignedDirect(14882)
                .reference("kraken:taxid|333760|NC_001526.4")
                .referenceTaxid(333760)
                .referenceKmerCount(35940)
                .alternateKmerCount(182529)
                .RName("adjusted_kraken_taxid_333760_NC_001526.4")
                .startPos(1)
                .endPos(7857)
                .numReads(5090)
                .covBases(6538)
                .coverage(83.2124)
                .meanDepth(90.4429)
                .meanBaseQ(35.9)
                .meanMapQ(59.7)
                .integrations(2)
                .qcStatus(VirusBreakendQCStatus.NO_ABNORMALITIES)
                .build();

        databaseAccess.writeVirusBreakend(TEST_SAMPLE_ID, List.of(virusBreakend));

        List<VirusbreakendRecord> virusBreakendRecords = fetchTable(Tables.VIRUSBREAKEND, VirusbreakendRecord.class);
        assertEquals(1, virusBreakendRecords.size());
    }
}
