package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.hartwig.hmftools.datamodel.purple.PurpleQC;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleQCFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleQCInterpretationTest
{
    @Test
    public void shouldInterpretContaminatedIfOneStateImpliesContamination()
    {
        assertTrue(PurpleQCInterpretation.isContaminated(create(PurpleQCStatus.FAIL_CONTAMINATION)));

        assertFalse(PurpleQCInterpretation.isContaminated(create(PurpleQCStatus.FAIL_NO_TUMOR)));

        assertTrue(PurpleQCInterpretation.isContaminated(create(PurpleQCStatus.PASS, PurpleQCStatus.FAIL_CONTAMINATION)));
    }

    @Test
    public void shouldInterpretAsFailIfAtLeastOneQCStateIsFail()
    {
        assertTrue(PurpleQCInterpretation.isQCFail(create(PurpleQCStatus.FAIL_NO_TUMOR)));
        assertTrue(PurpleQCInterpretation.isQCFail(create(PurpleQCStatus.FAIL_CONTAMINATION)));

        assertFalse(PurpleQCInterpretation.isQCFail(create(PurpleQCStatus.PASS)));

        assertTrue(PurpleQCInterpretation.isQCFail(create(PurpleQCStatus.PASS, PurpleQCStatus.FAIL_NO_TUMOR)));
    }

    @NotNull
    private static PurpleQC create(@NotNull PurpleQCStatus... states)
    {
        return TestPurpleQCFactory.builder().addAllStatus(Set.of(states)).build();
    }
}