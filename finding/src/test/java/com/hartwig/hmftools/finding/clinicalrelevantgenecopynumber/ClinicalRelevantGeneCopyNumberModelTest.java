package com.hartwig.hmftools.finding.clinicalrelevantgenecopynumber;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class ClinicalRelevantGeneCopyNumberModelTest
{
    @Test
    public void canExtractClinicalRelevantGeneCopyNumbers() {
        List<String> clinicalRelevantGeneCopyNumbersList = Lists.newArrayList();
        clinicalRelevantGeneCopyNumbersList.add("EGFR");
        ClinicalRelevantGeneCopyNumberModel clinicalRelevantCopyNumbersModel =
                new ClinicalRelevantGeneCopyNumberModel(clinicalRelevantGeneCopyNumbersList);

        assertTrue(clinicalRelevantCopyNumbersModel.findClinicalRelevantCopyNumberGene("EGFR"));
        assertFalse(clinicalRelevantCopyNumbersModel.findClinicalRelevantCopyNumberGene("BRAF"));
    }
}