package com.hartwig.hmftools.serve.vicc;

import static com.hartwig.hmftools.serve.vicc.ViccAnnotatedVCFChecker.extractGene;
import static com.hartwig.hmftools.serve.vicc.ViccAnnotatedVCFChecker.extractProteinAnnotation;
import static com.hartwig.hmftools.serve.vicc.ViccAnnotatedVCFChecker.extractTranscript;

import org.junit.Test;

public class ViccAnnotatedVCFCheckerTest {

    @Test
    public void testme() {
        String feature="PIK3CD:p.N334K - ENST00000377346";
        System.out.println(extractGene(feature));
        System.out.println(extractProteinAnnotation(feature));
        System.out.println(extractTranscript(feature));
    }

}