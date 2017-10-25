package com.hartwig.hmftools.patientdb.data;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class BiopsyTreatmentDataTest {

    @Test
    public void canGenerateCorrectTreatmentName() {
        final List<BiopsyTreatmentDrugData> drugs =
                Lists.newArrayList(drugWithName("drugB"), drugWithName("DrugC"), drugWithName("drugA"), drugWithName(null));

        final BiopsyTreatmentData data = withDrugs(drugs);

        assertEquals("DrugA/DrugB/DrugC/NULL", data.treatmentName());
    }

    @Test
    public void canGenerateCorrectTreatmentType() {
        final List<BiopsyTreatmentDrugData> noTypes = Lists.newArrayList(drugWithType(null));
        assertEquals("NULL", withDrugs(noTypes).type());

        final List<BiopsyTreatmentDrugData> simpleType = Lists.newArrayList(drugWithType("simple"), drugWithType("simple"));
        assertEquals("simple", withDrugs(simpleType).type());

        final List<BiopsyTreatmentDrugData> combiType = Lists.newArrayList(drugWithType("complex1"), drugWithType("complex2"));
        assertEquals("Combination", withDrugs(combiType).type());
    }

    @NotNull
    private static BiopsyTreatmentData withDrugs(@NotNull final List<BiopsyTreatmentDrugData> drugs) {
        return ImmutableBiopsyTreatmentData.of(1, null, null, null, drugs, null, Strings.EMPTY, Strings.EMPTY);
    }

    @NotNull
    private static BiopsyTreatmentDrugData drugWithName(@Nullable String name) {
        return ImmutableBiopsyTreatmentDrugData.of(name, null, null, null);
    }

    @NotNull
    private static BiopsyTreatmentDrugData drugWithType(@Nullable String type) {
        return ImmutableBiopsyTreatmentDrugData.of(null, type, null, null);
    }
}