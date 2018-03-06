package com.hartwig.hmftools.patientdb.data;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.time.LocalDate;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class TreatmentDataTest {

    @Test
    public void canGenerateCorrectTreatmentName() {
        List<DrugData> drugs =
                Lists.newArrayList(drugWithName("drugB"), drugWithName("DrugC"), drugWithName("drugA"), drugWithName(null));

        TreatmentData data = withDrugs(drugs);

        assertEquals("DrugA/DrugB/DrugC", data.treatmentName());
    }

    @Test
    public void canGenerateCorrectTreatmentType() {
        List<DrugData> noTypes = Lists.newArrayList(drugWithType(null));
        assertNull(withDrugs(noTypes).type());

        List<DrugData> simpleType = Lists.newArrayList(drugWithType("simple"), drugWithType("simple"));
        assertEquals("simple", withDrugs(simpleType).type());

        List<DrugData> combiType = Lists.newArrayList(drugWithType("complex1"), drugWithType("complex2"));
        assertEquals(TreatmentData.COMBI_THERAPY, withDrugs(combiType).type());
    }

    @Test
    public void sortsCorrectly() {
        TreatmentData treatment2015 = withDrug(drugWithStartDate(LocalDate.parse("2015-01-01")));
        TreatmentData treatment2014 = withDrug(drugWithStartDate(LocalDate.parse("2014-01-01")));
        TreatmentData treatmentNull = withDrug(drugWithStartDate(null));
        TreatmentData treatment2016 = withDrug(drugWithStartDate(LocalDate.parse("2016-01-01")));

        List<TreatmentData> treatments = Lists.newArrayList(treatment2015, treatment2014, treatmentNull, treatment2016);

        Collections.sort(treatments);

        assertEquals(treatment2014, treatments.get(0));
        assertEquals(treatment2015, treatments.get(1));
        assertEquals(treatment2016, treatments.get(2));
        assertEquals(treatmentNull, treatments.get(3));
    }

    @NotNull
    private static TreatmentData withDrug(@NotNull DrugData drug) {
        return ImmutableBiopsyTreatmentData.of(1, null, Lists.newArrayList(drug), null, FormStatusState.UNKNOWN, false);
    }

    @NotNull
    private static TreatmentData withDrugs(@NotNull List<DrugData> drugs) {
        return ImmutableBiopsyTreatmentData.of(1, null, drugs, null, FormStatusState.UNKNOWN, false);
    }

    @NotNull
    private static DrugData drugWithName(@Nullable String name) {
        final List<CuratedTreatment> curatedDrugs =
                name == null ? Lists.newArrayList() : Lists.newArrayList(ImmutableCuratedTreatment.of(name, "", ""));
        return ImmutableDrugData.of(name, null, null, null, curatedDrugs);
    }

    @NotNull
    private static DrugData drugWithType(@Nullable String type) {
        final List<CuratedTreatment> curatedDrugs =
                type == null ? Lists.newArrayList() : Lists.newArrayList(ImmutableCuratedTreatment.of("", type, ""));
        return ImmutableDrugData.of(null, null, null, null, curatedDrugs);
    }

    @NotNull
    private static DrugData drugWithStartDate(@Nullable LocalDate startDate) {
        return ImmutableDrugData.of(null, startDate, null, null, Lists.newArrayList());
    }
}