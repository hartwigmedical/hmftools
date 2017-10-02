package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.assertEquals;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class LimsJsonDataTest {

    @Test
    public void labProcedures() throws Exception {
        final String NA = "N/A";
        final String completeSopString = "PREP013V23-QC037V20-SEQ008V25";

        final LimsJsonData missingPrepNum = buildLimsEntry("PREPV23-QC037V20-SEQ008V25");
        assertEquals(NA, missingPrepNum.labProcedures());

        final LimsJsonData missingPrepVersion = buildLimsEntry("PREP013V-QC037V20-SEQ008V25");
        assertEquals(NA, missingPrepVersion.labProcedures());

        final LimsJsonData missingQcNum = buildLimsEntry("PREP013V23-QCV20-SEQ008V25");
        assertEquals(NA, missingQcNum.labProcedures());

        final LimsJsonData missingQcVersion = buildLimsEntry("PREP013V23-QC037V-SEQ008V25");
        assertEquals(NA, missingQcVersion.labProcedures());

        final LimsJsonData missingSeqNum = buildLimsEntry("PREP013V23-QC037V20-SEQV25");
        assertEquals(NA, missingSeqNum.labProcedures());

        final LimsJsonData missingSeqVersion = buildLimsEntry("PREP013V23-QC037V20-SEQ008V");
        assertEquals(NA, missingSeqVersion.labProcedures());

        final LimsJsonData completeEntry = buildLimsEntry(completeSopString);
        assertEquals(completeSopString, completeEntry.labProcedures());
    }

    @NotNull
    private static LimsJsonData buildLimsEntry(@NotNull final String labSopVersions) {
        return ImmutableLimsJsonData.builder()
                .sampleSource("")
                .patient("")
                .sampleName("")
                .sampleBarcode("")
                .arrivalDateString("")
                .samplingDateString("")
                .tumorPercentage("50")
                .labSopVersions(labSopVersions)
                .build();
    }
}