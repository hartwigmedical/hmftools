package com.hartwig.hmftools.patientreporter.copynumber;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.slicing.GenomeRegion;
import com.hartwig.hmftools.patientreporter.slicing.HMFSlicingAnnotation;
import com.hartwig.hmftools.patientreporter.slicing.HMFSlicingAnnotationFactory;

import org.junit.Test;

public class CopyNumberAnalyzerTest {

    private static final String CHROMOSOME = "X";
    private static final double EPSILON = 1.0e-10;

    @Test
    public void worksAsExpected() {
        final GenomeRegion first = new GenomeRegion(CHROMOSOME, 101, 200);
        final GenomeRegion second = new GenomeRegion(CHROMOSOME, 301, 400);
        final GenomeRegion third = new GenomeRegion(CHROMOSOME, 601, 700);

        final HMFSlicingAnnotation firstAnnotation = HMFSlicingAnnotationFactory.create("TRANS1", 1, "GENE1");
        final HMFSlicingAnnotation secondAnnotation = HMFSlicingAnnotationFactory.create("TRANS2", 1, "GENE2");
        final HMFSlicingAnnotation thirdAnnotation = HMFSlicingAnnotationFactory.create("TRANS3", 1, "GENE3");

        final Map<GenomeRegion, HMFSlicingAnnotation> annotations = Maps.newHashMap();
        annotations.put(first, firstAnnotation);
        annotations.put(second, secondAnnotation);
        annotations.put(third, thirdAnnotation);

        final List<CopyNumber> copyNumbers = Lists.newArrayList();
        copyNumbers.add(new CopyNumber(CHROMOSOME, 11, 20, 1));
        copyNumbers.add(new CopyNumber(CHROMOSOME, 121, 170, 4));
        copyNumbers.add(new CopyNumber(CHROMOSOME, 191, 260, 5));
        copyNumbers.add(new CopyNumber(CHROMOSOME, 291, 500, 0));

        final CopyNumberAnalyzer copyNumberAnalyzer = new CopyNumberAnalyzer(annotations);
        final CopyNumberAnalysis analysis = copyNumberAnalyzer.run(copyNumbers);

        final List<CopyNumberReport> findings = analysis.findings();
        assertEquals(1, findings.size());
        assertEquals(secondAnnotation.gene(), findings.get(0).gene());
        assertEquals(secondAnnotation.transcript(), findings.get(0).transcript());
        assertEquals(0, findings.get(0).copyNumber());

        final CopyNumberStats firstStat = analysis.stats().get(first);
        assertEquals(2, firstStat.min());
        assertEquals(5, firstStat.max());
        assertEquals(3.3, firstStat.mean(), EPSILON);

        final CopyNumberStats secondStat = analysis.stats().get(second);
        assertEquals(0, secondStat.min());
        assertEquals(0, secondStat.max());
        assertEquals(0D, secondStat.mean(), EPSILON);

        final CopyNumberStats thirdStat = analysis.stats().get(third);
        assertEquals(2, thirdStat.min());
        assertEquals(2, thirdStat.max());
        assertEquals(2D, thirdStat.mean(), EPSILON);
    }
}