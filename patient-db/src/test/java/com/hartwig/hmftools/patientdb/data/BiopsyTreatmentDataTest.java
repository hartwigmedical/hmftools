package com.hartwig.hmftools.patientdb.data;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class BiopsyTreatmentDataTest {

    @Test
    public void canGenerateCorrectTreatmentName() {
        final List<BiopsyTreatmentDrugData> drugs =
                Lists.newArrayList(drugWithName("drugB"), drugWithName("DrugC"), drugWithName("drugA"), drugWithName(null));

        final BiopsyTreatmentData data = withDrugs(drugs);

        assertEquals("DrugA/DrugB/DrugC", data.treatmentName());
    }

    @Test
    public void canGenerateCorrectTreatmentType() {
        final List<BiopsyTreatmentDrugData> noTypes = Lists.newArrayList(drugWithType(null));
        assertNull(withDrugs(noTypes).type());

        final List<BiopsyTreatmentDrugData> simpleType = Lists.newArrayList(drugWithType("simple"), drugWithType("simple"));
        assertEquals("simple", withDrugs(simpleType).type());

        final List<BiopsyTreatmentDrugData> combiType = Lists.newArrayList(drugWithType("complex1"), drugWithType("complex2"));
        assertEquals(BiopsyTreatmentData.COMBI_THERAPY, withDrugs(combiType).type());
    }

    @NotNull
    private static BiopsyTreatmentData withDrugs(@NotNull final List<BiopsyTreatmentDrugData> drugs) {
        return ImmutableBiopsyTreatmentData.of(1, null, null, null, drugs, null, FormStatusState.UNKNOWN, false);
    }

    @NotNull
    private static BiopsyTreatmentDrugData drugWithName(@Nullable String name) {
        final List<CuratedTreatment> curatedDrugs =
                name == null ? Lists.newArrayList() : Lists.newArrayList(ImmutableCuratedTreatment.of(name, "", ""));
        return ImmutableBiopsyTreatmentDrugData.of(name, null, null, curatedDrugs);
    }

    @NotNull
    private static BiopsyTreatmentDrugData drugWithType(@Nullable String type) {
        final List<CuratedTreatment> curatedDrugs =
                type == null ? Lists.newArrayList() : Lists.newArrayList(ImmutableCuratedTreatment.of("", type, ""));
        return ImmutableBiopsyTreatmentDrugData.of(null, null, null, curatedDrugs);
    }
}