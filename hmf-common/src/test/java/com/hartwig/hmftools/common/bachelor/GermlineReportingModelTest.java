package com.hartwig.hmftools.common.bachelor;

import static org.junit.Assert.*;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bachelor.GermlineReportingModel;

import org.junit.Test;

public class GermlineReportingModelTest {

    @Test
    public void returnFalseWhenGeneIsNotKnown() {
        GermlineReportingModel empty = new GermlineReportingModel(Maps.newHashMap());
        assertFalse(empty.notifyAboutGene("DoesNotExist"));
    }
}