package com.hartwig.hmftools.patientdb.data;

import static org.junit.Assert.assertEquals;

import java.time.LocalDate;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class BiopsyTreatmentResponseDataTest {

    @Test
    public void sortsCorrectly() {
        BiopsyTreatmentResponseData response2017 = responseOnDate(LocalDate.parse("2017-01-01"));
        BiopsyTreatmentResponseData responseNull = responseOnDate(null);
        BiopsyTreatmentResponseData response2015 = responseOnDate(LocalDate.parse("2015-01-01"));
        BiopsyTreatmentResponseData response2016 = responseOnDate(LocalDate.parse("2016-01-01"));

        List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(response2017, responseNull, response2015, response2016);

        Collections.sort(responses);
        assertEquals(response2015, responses.get(0));
        assertEquals(response2016, responses.get(1));
        assertEquals(response2017, responses.get(2));
        assertEquals(responseNull, responses.get(3));
    }

    @NotNull
    private static BiopsyTreatmentResponseData responseOnDate(@Nullable LocalDate date) {
        return ImmutableBiopsyTreatmentResponseData.builder()
                .formStatus(FormStatusState.UNKNOWN)
                .formLocked(false)
                .responseDate(date)
                .build();
    }
}