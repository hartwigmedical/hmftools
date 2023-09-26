package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.hartwig.hmftools.datamodel.purple.PurpleQC;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleQCFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleQualityInterpretationTest
{
    @Test
    public void shouldInterpretAsFailIfAtLeastOneQCStateIsFail()
    {
        assertTrue(PurpleQualityInterpretation.isQCFail(create(PurpleQCStatus.FAIL_NO_TUMOR)));
        assertTrue(PurpleQualityInterpretation.isQCFail(create(PurpleQCStatus.FAIL_CONTAMINATION)));

        assertFalse(PurpleQualityInterpretation.isQCFail(create(PurpleQCStatus.PASS)));

        assertTrue(PurpleQualityInterpretation.isQCFail(create(PurpleQCStatus.PASS, PurpleQCStatus.FAIL_NO_TUMOR)));
    }

    @NotNull
    private static PurpleQC create(@NotNull PurpleQCStatus... states)
    {
        return TestPurpleQCFactory.builder().addAllStatus(Set.of(states)).build();
    }
}