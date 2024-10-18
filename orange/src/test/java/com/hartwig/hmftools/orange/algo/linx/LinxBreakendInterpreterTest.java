package com.hartwig.hmftools.orange.algo.linx;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.gene.TranscriptRegionType;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.orange.algo.pave.TestEnsemblDataCacheFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class LinxBreakendInterpreterTest
{

    private static final String GENE = "gene";
    private static final String VCF_ID = "vcfId";
    private static final double EPSILON = 0.000001D;

    @Test
    public void canInterpretValidData()
    {
        List<StructuralVariant> structuralVariants = List.of(createStructuralVariant());

        List<LinxSvAnnotation> linxSvAnnotations = List.of(createSvAnnotation());
        List<LinxBreakend> breakends = createBreakends();

        LinxBreakendInterpreter interpreter = createInterpreter(structuralVariants, linxSvAnnotations);
        com.hartwig.hmftools.datamodel.linx.LinxBreakend left = interpreter.interpret(breakends.get(0));
        assertEquals("7", left.chromosome());
        assertEquals("q34", left.chromosomeBand());
        assertEquals(LinxBreakendType.DUP, left.type());
        assertEquals(-1, left.orientation());
        assertEquals(1.5D, left.junctionCopyNumber(), EPSILON);

        com.hartwig.hmftools.datamodel.linx.LinxBreakend right = interpreter.interpret(breakends.get(1));
        assertEquals("7", right.chromosome());
        assertEquals("q34", right.chromosomeBand());
        assertEquals(LinxBreakendType.DUP, right.type());
        assertEquals(1, right.orientation());
        assertEquals(1.5D, right.junctionCopyNumber(), EPSILON);
    }

    @Test
    public void canFallbackWhenNoSvAnnotation()
    {
        List<StructuralVariant> structuralVariants = List.of(createStructuralVariant());

        List<LinxSvAnnotation> linxSvAnnotations = List.of();
        List<LinxBreakend> breakends = createBreakends();

        LinxBreakendInterpreter interpreter = createInterpreter(structuralVariants, linxSvAnnotations);
        com.hartwig.hmftools.datamodel.linx.LinxBreakend left = interpreter.interpret(breakends.get(0));
        assertEquals("", left.chromosome());
        assertEquals("q34", left.chromosomeBand());
        assertEquals(LinxBreakendType.BND, left.type());
        assertEquals(0, left.orientation());
        assertEquals(0D, left.junctionCopyNumber(), EPSILON);

        com.hartwig.hmftools.datamodel.linx.LinxBreakend right = interpreter.interpret(breakends.get(1));
        assertEquals("", right.chromosome());
        assertEquals("q34", right.chromosomeBand());
        assertEquals(LinxBreakendType.BND, right.type());
        assertEquals(0, right.orientation());
        assertEquals(0D, right.junctionCopyNumber(), EPSILON);
    }

    @Test
    public void canFallbackWhenNoStructuralVariant()
    {
        List<StructuralVariant> structuralVariants = List.of();

        List<LinxSvAnnotation> linxSvAnnotations = List.of(createSvAnnotation());
        List<LinxBreakend> breakends = createBreakends();

        LinxBreakendInterpreter interpreter = createInterpreter(structuralVariants, linxSvAnnotations);
        com.hartwig.hmftools.datamodel.linx.LinxBreakend left = interpreter.interpret(breakends.get(0));
        assertEquals("", left.chromosome());
        assertEquals("q34", left.chromosomeBand());
        assertEquals(LinxBreakendType.BND, left.type());
        assertEquals(0, left.orientation());
        assertEquals(1.5D, left.junctionCopyNumber(), 0.0D);

        com.hartwig.hmftools.datamodel.linx.LinxBreakend right = interpreter.interpret(breakends.get(1));
        assertEquals("", right.chromosome());
        assertEquals("q34", right.chromosomeBand());
        assertEquals(LinxBreakendType.BND, right.type());
        assertEquals(0, right.orientation());
        assertEquals(1.5D, right.junctionCopyNumber(), 0.0D);
    }

    @Test
    public void canComputeJunctionCopyNumber()
    {
        LinxSvAnnotation svAnnotation = createSvAnnotation();

        double result = LinxBreakendInterpreter.junctionCopyNumber(svAnnotation);
        assertEquals(1.5, result, EPSILON);
    }

    @NotNull
    private static StructuralVariant createStructuralVariant()
    {
        return PurpleTestUtils.createStructuralVariant("7", 50, "7", 200, StructuralVariantType.DUP, 1.0, 1.0)
                .id(VCF_ID)
                .build();
    }

    @NotNull
    private static List<LinxBreakend> createBreakends()
    {
        return List.of(
                createBreakend(1, true),
                createBreakend(2, false)
        );
    }

    @NotNull
    private static LinxBreakend createBreakend(int id, boolean isStart)
    {
        return LinxTestFactory.breakendBuilder()
                .id(id)
                .reportedDisruption(true)
                .gene(GENE)
                .isStart(isStart)
                .svId(1)
                .regionType(TranscriptRegionType.INTRONIC)
                .build();
    }

    @NotNull
    private static LinxSvAnnotation createSvAnnotation()
    {
        return LinxTestFactory.svAnnotationBuilder()
                .svId(1)
                .vcfId(VCF_ID)
                .junctionCopyNumberMin(1.0D)
                .junctionCopyNumberMax(2.0D)
                .build();
    }

    @NotNull
    private static LinxBreakendInterpreter createInterpreter(@NotNull List<StructuralVariant> structuralVariants,
            @NotNull List<LinxSvAnnotation> linxSvAnnotations)
    {
        return new LinxBreakendInterpreter(structuralVariants, linxSvAnnotations, TestEnsemblDataCacheFactory.loadTestCache());
    }
}