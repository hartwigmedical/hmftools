package com.hartwig.hmftools.patientdb.clinical.datamodel;

import static com.hartwig.hmftools.patientdb.clinical.datamodel.ClinicalTestFactory.biopsyTreatmentResponseBuilder;

import static org.junit.Assert.assertEquals;

import java.time.LocalDate;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class BiopsyTreatmentResponseDataTest {

    @Test
    public void sortsCorrectly() {
        BiopsyTreatmentResponseData response2017 = biopsyTreatmentResponseBuilder().responseDate(LocalDate.parse("2017-01-01")).build();
        BiopsyTreatmentResponseData responseNull = biopsyTreatmentResponseBuilder().responseDate(null).build();
        BiopsyTreatmentResponseData response2015 = biopsyTreatmentResponseBuilder().responseDate(LocalDate.parse("2015-01-01")).build();
        BiopsyTreatmentResponseData response2016 = biopsyTreatmentResponseBuilder().responseDate(LocalDate.parse("2016-01-01")).build();

        List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(response2017, responseNull, response2015, response2016);

        Collections.sort(responses);
        assertEquals(response2015, responses.get(0));
        assertEquals(response2016, responses.get(1));
        assertEquals(response2017, responses.get(2));
        assertEquals(responseNull, responses.get(3));
    }
}