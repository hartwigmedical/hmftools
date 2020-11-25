package com.hartwig.hmftools.common.clinvar;

import static com.hartwig.hmftools.common.clinvar.ClinvarPathogenicity.BENIGN;
import static com.hartwig.hmftools.common.clinvar.ClinvarPathogenicity.CLINVAR_STR_BENIGN;
import static com.hartwig.hmftools.common.clinvar.ClinvarPathogenicity.CLINVAR_STR_CONFLICTING;
import static com.hartwig.hmftools.common.clinvar.ClinvarPathogenicity.CLINVAR_STR_LIKELY_BENIGN;
import static com.hartwig.hmftools.common.clinvar.ClinvarPathogenicity.CLINVAR_STR_LIKELY_PATHOGENIC;
import static com.hartwig.hmftools.common.clinvar.ClinvarPathogenicity.CLINVAR_STR_PATHOGENIC;
import static com.hartwig.hmftools.common.clinvar.ClinvarPathogenicity.CONFLICTING;
import static com.hartwig.hmftools.common.clinvar.ClinvarPathogenicity.LIKELY_BENIGN;
import static com.hartwig.hmftools.common.clinvar.ClinvarPathogenicity.LIKELY_PATHOGENIC;
import static com.hartwig.hmftools.common.clinvar.ClinvarPathogenicity.PATHOGENIC;
import static com.hartwig.hmftools.common.clinvar.ClinvarPathogenicity.UNKNOWN;
import static com.hartwig.hmftools.common.clinvar.ClinvarPathogenicity.fromClinvarAnnotation;

import static org.junit.Assert.assertEquals;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class ClinvarPathogenicityTest {

    @Test
    public void testFromAnnotation() {
        assertEquals(UNKNOWN, fromClinvarAnnotation(Strings.EMPTY, Strings.EMPTY));

        assertEquals(PATHOGENIC, fromClinvarAnnotation(CLINVAR_STR_PATHOGENIC, Strings.EMPTY));

        assertEquals(PATHOGENIC, fromClinvarAnnotation(CLINVAR_STR_PATHOGENIC, Strings.EMPTY));
        assertEquals(PATHOGENIC, fromClinvarAnnotation(CLINVAR_STR_PATHOGENIC + CLINVAR_STR_LIKELY_PATHOGENIC, Strings.EMPTY));
        assertEquals(PATHOGENIC, fromClinvarAnnotation(CLINVAR_STR_CONFLICTING, CLINVAR_STR_PATHOGENIC + CLINVAR_STR_LIKELY_PATHOGENIC));

        assertEquals(LIKELY_PATHOGENIC, fromClinvarAnnotation(CLINVAR_STR_LIKELY_PATHOGENIC, Strings.EMPTY));
        assertEquals(LIKELY_PATHOGENIC, fromClinvarAnnotation(CLINVAR_STR_CONFLICTING,  CLINVAR_STR_LIKELY_PATHOGENIC));

        assertEquals(BENIGN, fromClinvarAnnotation(CLINVAR_STR_BENIGN, Strings.EMPTY));
        assertEquals(BENIGN, fromClinvarAnnotation(CLINVAR_STR_BENIGN + CLINVAR_STR_LIKELY_BENIGN, Strings.EMPTY));
        assertEquals(BENIGN, fromClinvarAnnotation(CLINVAR_STR_CONFLICTING, CLINVAR_STR_BENIGN + CLINVAR_STR_LIKELY_BENIGN));

        assertEquals(LIKELY_BENIGN, fromClinvarAnnotation(CLINVAR_STR_LIKELY_BENIGN, Strings.EMPTY));
        assertEquals(LIKELY_BENIGN, fromClinvarAnnotation(CLINVAR_STR_CONFLICTING,  CLINVAR_STR_LIKELY_BENIGN));

        assertEquals(CONFLICTING, fromClinvarAnnotation(CLINVAR_STR_CONFLICTING,  CLINVAR_STR_PATHOGENIC + CLINVAR_STR_BENIGN));
        assertEquals(CONFLICTING, fromClinvarAnnotation(CLINVAR_STR_CONFLICTING,  CLINVAR_STR_PATHOGENIC + CLINVAR_STR_LIKELY_BENIGN));
        assertEquals(CONFLICTING, fromClinvarAnnotation(CLINVAR_STR_CONFLICTING,  CLINVAR_STR_LIKELY_PATHOGENIC + CLINVAR_STR_BENIGN));
        assertEquals(CONFLICTING, fromClinvarAnnotation(CLINVAR_STR_CONFLICTING,  CLINVAR_STR_LIKELY_PATHOGENIC + CLINVAR_STR_LIKELY_BENIGN));

    }

}
