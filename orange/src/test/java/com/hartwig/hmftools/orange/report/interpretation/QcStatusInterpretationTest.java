package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.hartwig.hmftools.datamodel.purple.PurpleQC;
import com.hartwig.hmftools.orange.algo.QcStatusInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleQCFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class QcStatusInterpretationTest
{
    @Test
    public void shouldInterpretAsFailIfAtLeastOneQCStateIsFail()
    {
        assertTrue(QcStatusInterpretation.hasPurpleFail(create(PurpleQCStatus.FAIL_NO_TUMOR)));
        assertTrue(QcStatusInterpretation.hasPurpleFail(create(PurpleQCStatus.FAIL_CONTAMINATION)));

        assertFalse(QcStatusInterpretation.hasPurpleFail(create(PurpleQCStatus.PASS)));

        assertTrue(QcStatusInterpretation.hasPurpleFail(create(PurpleQCStatus.PASS, PurpleQCStatus.FAIL_NO_TUMOR)));
    }

    @Test
    public void shouldInterpretContaminatedIfOneStateImpliesContamination()
    {
        assertTrue(QcStatusInterpretation.hasTumorContaminated(create(PurpleQCStatus.FAIL_CONTAMINATION)));

        assertFalse(QcStatusInterpretation.hasTumorContaminated(create(PurpleQCStatus.FAIL_NO_TUMOR)));

        assertTrue(QcStatusInterpretation.hasTumorContaminated(create(PurpleQCStatus.PASS, PurpleQCStatus.FAIL_CONTAMINATION)));
    }

    @Test
    public void shouldInterpretFailNoTumorIfOneStateImpliesFailNoTumor()
    {
        assertTrue(QcStatusInterpretation.isFailNoTumor(create(PurpleQCStatus.FAIL_NO_TUMOR)));

        assertFalse(QcStatusInterpretation.isFailNoTumor(create(PurpleQCStatus.FAIL_CONTAMINATION)));

        assertTrue(QcStatusInterpretation.isFailNoTumor(create(PurpleQCStatus.PASS, PurpleQCStatus.FAIL_NO_TUMOR)));
    }

    @NotNull
    private static PurpleQC create(@NotNull PurpleQCStatus... states)
    {
        return TestPurpleQCFactory.builder().addAllStatus(Set.of(states)).build();
    }
}