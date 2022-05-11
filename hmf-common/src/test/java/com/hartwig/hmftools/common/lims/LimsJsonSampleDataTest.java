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

        String completeSopStringOld = "PREP013V23-QC037V20-SEQ008V25";
        LimsJsonSampleData completeEntryOld = buildLimsSampleEntry(completeSopStringOld);
        assertEquals(completeSopStringOld, completeEntryOld.labProcedures());

        String completeSopStringNew = "PREP013V23-ENR11V12-QC037V20-SEQ008V25";
        LimsJsonSampleData completeEntryNew = buildLimsSampleEntry(completeSopStringNew);
        assertEquals(completeSopStringNew, completeEntryNew.labProcedures());

        String completeSopStringNewWithoutENR = "PREP013V23-ENRnaVna-QC037V20-SEQ008V25";
        String completeSopStringUsable = "PREP013V23-QC037V20-SEQ008V25";
        LimsJsonSampleData completeSopStringNewWithoutENRCode = buildLimsSampleEntry(completeSopStringNewWithoutENR);
        assertEquals(completeSopStringUsable, completeSopStringNewWithoutENRCode.labProcedures());
    }

    @NotNull
    private static LimsJsonSampleData buildLimsSampleEntry(@NotNull String labSopVersions) {
        return createLimsSampleDataBuilder().labSopVersions(labSopVersions).build();
    }
}