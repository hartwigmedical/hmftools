package com.hartwig.hmftools.common.lims;

import static com.hartwig.hmftools.common.lims.LimsTestUtil.createLimsSampleDataBuilder;

import static org.junit.Assert.assertEquals;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class LimsJsonSampleDataTest {

    @Test
    public void canConvertSOPsToLabProceduresString() {
        LimsJsonSampleData missingPrepNum = buildLimsSampleEntry("PREPV23-QC037V20-SEQ008V25");
        assertEquals(Lims.NOT_AVAILABLE_STRING, missingPrepNum.labProcedures());

        LimsJsonSampleData missingPrepVersion = buildLimsSampleEntry("PREP013V-QC037V20-SEQ008V25");
        assertEquals(Lims.NOT_AVAILABLE_STRING, missingPrepVersion.labProcedures());

        LimsJsonSampleData missingQcNum = buildLimsSampleEntry("PREP013V23-QCV20-SEQ008V25");
        assertEquals(Lims.NOT_AVAILABLE_STRING, missingQcNum.labProcedures());

        LimsJsonSampleData missingQcVersion = buildLimsSampleEntry("PREP013V23-QC037V-SEQ008V25");
        assertEquals(Lims.NOT_AVAILABLE_STRING, missingQcVersion.labProcedures());

        LimsJsonSampleData missingSeqNum = buildLimsSampleEntry("PREP013V23-QC037V20-SEQV25");
        assertEquals(Lims.NOT_AVAILABLE_STRING, missingSeqNum.labProcedures());

        LimsJsonSampleData missingSeqVersion = buildLimsSampleEntry("PREP013V23-QC037V20-SEQ008V");
        assertEquals(Lims.NOT_AVAILABLE_STRING, missingSeqVersion.labProcedures());

        String completeSopString = "PREP013V23-QC037V20-SEQ008V25";
        LimsJsonSampleData completeEntry = buildLimsSampleEntry(completeSopString);
        assertEquals(completeSopString, completeEntry.labProcedures());
    }

    @NotNull
    private static LimsJsonSampleData buildLimsSampleEntry(@NotNull String labSopVersions) {
        return createLimsSampleDataBuilder().labSopVersions(labSopVersions).build();
    }
}