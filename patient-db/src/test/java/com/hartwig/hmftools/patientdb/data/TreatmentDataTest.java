package com.hartwig.hmftools.patientdb.data;

import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.biopsyTreatmentBuilder;
import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.drugBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.time.LocalDate;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class TreatmentDataTest {

    @Test
    public void canGenerateCorrectTreatmentName() {
        List<DrugData> drugs = Lists.newArrayList(drugWithName("DrugB"), drugWithName("DrugC"), drugWithName("drugA"), drugWithName(null));

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
        return biopsyTreatmentBuilder().addDrugs(drug).build();
    }

    @NotNull
    private static TreatmentData withDrugs(@NotNull List<DrugData> drugs) {
        return biopsyTreatmentBuilder().addAllDrugs(drugs).build();
    }

    @NotNull
    private static DrugData drugWithName(@Nullable String name) {
        return name != null
                ? drugBuilder().addCuratedTreatments(ImmutableCuratedTreatment.of(name, "", "")).build()
                : drugBuilder().build();
    }

    @NotNull
    private static DrugData drugWithType(@Nullable String type) {
        return type != null
                ? drugBuilder().addCuratedTreatments(ImmutableCuratedTreatment.of("", type, "")).build()
                : drugBuilder().build();
    }

    @NotNull
    private static DrugData drugWithStartDate(@Nullable LocalDate startDate) {
        return drugBuilder().startDate(startDate).build();
    }
}