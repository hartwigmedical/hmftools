package com.hartwig.hmftools.patientdb;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.ImmutableAnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.VirusannotationRecord;

import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class LoadVirusInterpreterTest extends DatabaseTestBase
{
    @Test
    public void canWriteVirusInterpreterData()
    {
        AnnotatedVirus virus = ImmutableAnnotatedVirus.builder()
                .taxid(333760)
                .name("Human papillomavirus 16")
                .qcStatus(VirusBreakendQCStatus.NO_ABNORMALITIES)
                .integrations(5)
                .interpretation(VirusType.HPV)
                .percentageCovered(96.4963)
                .meanCoverage(270.166)
                .expectedClonalCoverage(19.54551099)
                .reported(true)
                .blacklisted(false)
                .virusDriverLikelihoodType(VirusLikelihoodType.HIGH)
                .build();

        databaseAccess.writeVirusInterpreter(TEST_SAMPLE_ID, List.of(virus));

        List<VirusannotationRecord> records = fetchTable(Tables.VIRUSANNOTATION, VirusannotationRecord.class);
        assertEquals(1, records.size());
    }
}
