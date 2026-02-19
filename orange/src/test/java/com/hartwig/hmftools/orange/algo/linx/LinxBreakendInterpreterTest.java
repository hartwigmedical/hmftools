package com.hartwig.hmftools.orange.algo.linx;

import static com.hartwig.hmftools.orange.TestDataUtils.CYTO_BANDS;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.gene.TranscriptRegionType;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.orange.TestDataUtils;

import org.junit.Test;

public class LinxBreakendInterpreterTest
{
    private static final String GENE = "gene";
    private static final String VCF_ID = "vcfId";
    private static final double EPSILON = 0.000001D;

    @Test
    public void canInterpretValidData()
    {
        List<LinxSvAnnotation> linxSvAnnotations = List.of(createSvAnnotation());
        List<LinxBreakend> breakends = createBreakends();

        LinxBreakendInterpreter interpreter = createInterpreter(linxSvAnnotations);
        com.hartwig.hmftools.datamodel.linx.LinxBreakend left = interpreter.interpret(breakends.get(0));
        assertEquals("1", left.chromosome());
        assertEquals("p36.33", left.chromosomeBand());
        assertEquals(LinxBreakendType.DEL, left.type());
        assertEquals(1, left.orientation());
        assertEquals(1.5D, left.junctionCopyNumber(), EPSILON);

        com.hartwig.hmftools.datamodel.linx.LinxBreakend right = interpreter.interpret(breakends.get(1));
        assertEquals("1", right.chromosome());
        assertEquals("p36.33", right.chromosomeBand());
        assertEquals(LinxBreakendType.DEL, right.type());
        assertEquals(-1, right.orientation());
        assertEquals(1.5D, right.junctionCopyNumber(), EPSILON);
    }

    @Test
    public void canComputeJunctionCopyNumber()
    {
        LinxSvAnnotation svAnnotation = createSvAnnotation();

        double result = LinxBreakendInterpreter.junctionCopyNumber(svAnnotation);
        assertEquals(1.5, result, EPSILON);
    }

    private static List<LinxBreakend> createBreakends()
    {
        return List.of(
                createBreakend(1, true),
                createBreakend(2, false)
        );
    }

    private static LinxBreakend createBreakend(int id, boolean isStart)
    {
        return LinxTestFactory.breakendBuilder()
                .id(id)
                .reportedStatus(ReportedStatus.REPORTED)
                .gene(GENE)
                .isStart(isStart)
                .svId(1)
                .regionType(TranscriptRegionType.INTRONIC)
                .build();
    }

    private static LinxSvAnnotation createSvAnnotation()
    {
        return LinxTestFactory.svAnnotationBuilder()
                .svId(1)
                .vcfIdStart(VCF_ID)
                .vcfIdEnd(VCF_ID)
                .junctionCopyNumberMin(1.0D)
                .junctionCopyNumberMax(2.0D)
                .build();
    }

    private static LinxBreakendInterpreter createInterpreter(final List<LinxSvAnnotation> linxSvAnnotations)
    {
        return new LinxBreakendInterpreter(linxSvAnnotations, CYTO_BANDS);
    }
}