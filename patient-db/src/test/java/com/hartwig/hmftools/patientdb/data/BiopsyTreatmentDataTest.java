package com.hartwig.hmftools.patientdb.data;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.time.LocalDate;
import java.util.List;
import java.util.SortedSet;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class BiopsyTreatmentDataTest {

    @Test
    public void canGenerateCorrectTreatmentName() {
        List<BiopsyTreatmentDrugData> drugs =
                Lists.newArrayList(drugWithName("drugB"), drugWithName("DrugC"), drugWithName("drugA"), drugWithName(null));

        BiopsyTreatmentData data = withDrugs(drugs);

        assertEquals("DrugA/DrugB/DrugC", data.treatmentName());
    }

    @Test
    public void canGenerateCorrectTreatmentType() {
        List<BiopsyTreatmentDrugData> noTypes = Lists.newArrayList(drugWithType(null));
        assertNull(withDrugs(noTypes).type());

        List<BiopsyTreatmentDrugData> simpleType = Lists.newArrayList(drugWithType("simple"), drugWithType("simple"));
        assertEquals("simple", withDrugs(simpleType).type());

        List<BiopsyTreatmentDrugData> combiType = Lists.newArrayList(drugWithType("complex1"), drugWithType("complex2"));
        assertEquals(BiopsyTreatmentData.COMBI_THERAPY, withDrugs(combiType).type());
    }

    @Test
    public void sortsCorrectly() {
        BiopsyTreatmentData treatment2015 = withDrug(drugWithStartDate(LocalDate.parse("2015-01-01")));
        BiopsyTreatmentData treatment2014 = withDrug(drugWithStartDate(LocalDate.parse("2014-01-01")));
        BiopsyTreatmentData treatmentNull = withDrug(drugWithStartDate(null));
        BiopsyTreatmentData treatment2016 = withDrug(drugWithStartDate(LocalDate.parse("2016-01-01")));

        List<BiopsyTreatmentData> unsorted = Lists.newArrayList(treatment2015, treatment2014, treatmentNull, treatment2016);

        SortedSet<BiopsyTreatmentData> sorted = Sets.newTreeSet(unsorted);

        assertEquals(treatment2014, sorted.first());
        assertEquals(treatmentNull, sorted.last());
    }

    @NotNull
    private static BiopsyTreatmentData withDrug(@NotNull BiopsyTreatmentDrugData drug) {
        return ImmutableBiopsyTreatmentData.of(1, null, Lists.newArrayList(drug), null, FormStatusState.UNKNOWN, false);
    }

    @NotNull
    private static BiopsyTreatmentData withDrugs(@NotNull List<BiopsyTreatmentDrugData> drugs) {
        return ImmutableBiopsyTreatmentData.of(1, null, drugs, null, FormStatusState.UNKNOWN, false);
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

    @NotNull
    private static BiopsyTreatmentDrugData drugWithStartDate(@Nullable LocalDate startDate) {
        return ImmutableBiopsyTreatmentDrugData.of(null, startDate, null, Lists.newArrayList());
    }
}