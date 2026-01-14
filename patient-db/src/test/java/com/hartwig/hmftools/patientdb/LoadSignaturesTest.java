package com.hartwig.hmftools.patientdb;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.sigs.ImmutableSignatureAllocation;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.SignatureRecord;

import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class LoadSignaturesTest extends DatabaseTestBase
{
    @Test
    public void canWriteSignatures()
    {
        SignatureAllocation signatureAllocation = ImmutableSignatureAllocation.builder()
                .signature("Sig2")
                .allocation(3882.679)
                .percent(0.10253)
                .build();

        databaseAccess.writeSignatures(TEST_SAMPLE_ID, List.of(signatureAllocation));

        List<SignatureRecord> signatureRecords = fetchTable(Tables.SIGNATURE, SignatureRecord.class);
        assertEquals(1, signatureRecords.size());
    }
}
