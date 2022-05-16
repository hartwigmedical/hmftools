package com.hartwig.hmftools.patientdb.clinical.datamodel;

import static com.hartwig.hmftools.patientdb.clinical.datamodel.DatamodelTestFactory.biopsyTreatmentBuilder;
import static com.hartwig.hmftools.patientdb.clinical.datamodel.DatamodelTestFactory.drugBuilder;

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
        List<DrugData> drugs = Lists.newArrayList(drugWithName("DrugB"), drugWithName("DrugC"), drugWithName("DrugA"), drugWithName(null));

        TreatmentData data = withDrugs(drugs);

        assertEquals("DrugB/DrugC/DrugA", data.treatmentName());
    }

    @Test
    public void canGenerateCorrectTreatmentType() {
        List<DrugData> noTypes = Lists.newArrayList(drugWithType(null));
        assertNull(withDrugs(noTypes).consolidatedType());
        assertNull(withDrugs(noTypes).concatenatedType());

        List<DrugData> simpleType = Lists.newArrayList(drugWithType("simple"), drugWithType("simple"));
        assertEquals("simple", withDrugs(simpleType).consolidatedType());
        assertEquals("simple/simple", withDrugs(simpleType).concatenatedType());

        List<DrugData> combiType = Lists.newArrayList(drugWithType("complex1"), drugWithType("complex2"));
        assertEquals(TreatmentData.COMBI_TYPE, withDrugs(combiType).consolidatedType());
        assertEquals("complex1/complex2", withDrugs(combiType).concatenatedType());
    }

    @Test
    public void canGenerateCorrectTreatmentMechanism() {
        List<DrugData> noMechanism = Lists.newArrayList(drugWithMechanism(null));
        assertNull(withDrugs(noMechanism).consolidatedMechanism());
        assertNull(withDrugs(noMechanism).concatenatedMechanism());

        List<DrugData> simpleMechanism = Lists.newArrayList(drugWithMechanism("mechanism"), drugWithMechanism("mechanism"));
        assertEquals("mechanism", withDrugs(simpleMechanism).consolidatedMechanism());
        assertEquals("mechanism/mechanism", withDrugs(simpleMechanism).concatenatedMechanism());

        List<DrugData> combiMechanism = Lists.newArrayList(drugWithMechanism("mechanism1"), drugWithMechanism("mechanism2"));
        assertEquals(TreatmentData.COMBI_MECHANISM, withDrugs(combiMechanism).consolidatedMechanism());
        assertEquals("mechanism1/mechanism2", withDrugs(combiMechanism).concatenatedMechanism());
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
        return name != null ? drugBuilder().addCuratedDrugs(ImmutableCuratedDrug.of(name, "", "", "")).build()
                : drugBuilder().build();
    }

    @NotNull
    private static DrugData drugWithType(@Nullable String type) {
        return type != null ? drugBuilder().addCuratedDrugs(ImmutableCuratedDrug.of("", type, "", "")).build()
                : drugBuilder().build();
    }

    @NotNull
    private static DrugData drugWithMechanism(@Nullable String mechanism) {
        return mechanism != null
                ? drugBuilder().addCuratedDrugs(ImmutableCuratedDrug.of("", "", mechanism, "")).build()
                : drugBuilder().build();
    }

    @NotNull
    private static DrugData drugWithStartDate(@Nullable LocalDate startDate) {
        return drugBuilder().startDate(startDate).build();
    }
}