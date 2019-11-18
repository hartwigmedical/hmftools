package com.hartwig.hmftools.patientreporter.variants.germline;

import static org.junit.Assert.*;

import com.google.common.collect.Maps;

import org.junit.Test;

public class GermlineReportingModelTest {

    @Test
    public void returnFalseWhenGeneIsNotKnown() {
        GermlineReportingModel empty = new GermlineReportingModel(Maps.newHashMap());
        assertFalse(empty.notifyAboutGene("DoesNotExist"));
    }
}