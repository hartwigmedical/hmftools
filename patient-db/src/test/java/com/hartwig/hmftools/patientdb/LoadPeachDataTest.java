package com.hartwig.hmftools.patientdb;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.peach.ImmutablePeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.PeachgenotypeRecord;

import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class LoadPeachDataTest extends DatabaseTestBase
{
    @Test
    public void canWritePeachGenotype()
    {
        PeachGenotype genotype = ImmutablePeachGenotype.builder()
                .gene("UGT1A1")
                .allele("*1")
                .alleleCount(2)
                .function("Normal Function")
                .linkedDrugs("Irinotecan")
                .urlPrescriptionInfo("https://www.pharmgkb.org/guidelineAnnotation/PA166104951")
                .build();

        databaseAccess.writePeach(TEST_SAMPLE_ID, List.of(genotype));

        List<PeachgenotypeRecord> records = fetchTable(Tables.PEACHGENOTYPE, PeachgenotypeRecord.class);
        assertEquals(1, records.size());
    }
}
