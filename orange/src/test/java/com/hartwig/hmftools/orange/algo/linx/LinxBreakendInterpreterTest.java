package com.hartwig.hmftools.orange.algo.linx;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.orange.algo.pave.TestEnsemblDataCacheFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class LinxBreakendInterpreterTest
{

    private static final double EPSILON = 0.000001D;

    @Test
    public void canInterpretValidData()
    {
        List<StructuralVariant> structuralVariants = List.of(
                createStructuralVariant(StructuralVariantType.DUP)
        );

        List<LinxSvAnnotation> linxSvAnnotations = List.of(svAnnotation());
        List<LinxBreakend> breakends = createBreakends();

        LinxBreakendInterpreter interpreter = createInterpreter(structuralVariants, linxSvAnnotations);
        com.hartwig.hmftools.datamodel.linx.LinxBreakend left = interpreter.interpret(breakends.get(0));
        assertNotNull(left);
        assertEquals("7", left.chromosome());
        assertEquals("q34", left.chromosomeBand());
        assertEquals(LinxBreakendType.DUP, left.type());
        assertEquals(-1, left.orientation());
        assertEquals(1.5D, left.junctionCopyNumber(), EPSILON);

        com.hartwig.hmftools.datamodel.linx.LinxBreakend right = interpreter.interpret(breakends.get(1));
        assertNotNull(right);
        assertEquals("7", right.chromosome());
        assertEquals("q34", right.chromosomeBand());
        assertEquals(LinxBreakendType.DUP, right.type());
        assertEquals(1, right.orientation());
        assertEquals(1.5D, right.junctionCopyNumber(), EPSILON);
    }

    @Test
    public void fallbackWhenNoSvAnnotation()
    {
        List<StructuralVariant> structuralVariants = List.of(
                createStructuralVariant(StructuralVariantType.DUP)
        );

        List<LinxSvAnnotation> linxSvAnnotations = List.of();
        List<LinxBreakend> breakends = createBreakends();

        LinxBreakendInterpreter interpreter = createInterpreter(structuralVariants, linxSvAnnotations);
        com.hartwig.hmftools.datamodel.linx.LinxBreakend left = interpreter.interpret(breakends.get(0));
        assertFallbackBreakendAttributes(left);

        com.hartwig.hmftools.datamodel.linx.LinxBreakend right = interpreter.interpret(breakends.get(1));
        assertFallbackBreakendAttributes(right);
    }

    @Test
    public void fallbackWhenNoStructuralVariant()
    {
        List<StructuralVariant> structuralVariants = List.of();

        List<LinxSvAnnotation> linxSvAnnotations = List.of(svAnnotation());
        List<LinxBreakend> breakends = createBreakends();

        LinxBreakendInterpreter interpreter = createInterpreter(structuralVariants, linxSvAnnotations);
        com.hartwig.hmftools.datamodel.linx.LinxBreakend left = interpreter.interpret(breakends.get(0));
        assertFallbackBreakendAttributes(left);

        com.hartwig.hmftools.datamodel.linx.LinxBreakend right = interpreter.interpret(breakends.get(1));
        assertFallbackBreakendAttributes(right);
    }

    @Test
    public void fallbackWhenNoEndLeg()
    {
        List<StructuralVariant> structuralVariants = List.of(
                PurpleTestUtils.createStructuralVariantSingleBreakend("7", 50, 1.0)
                        .id("vcfId")
                        .build()
        );

        List<LinxSvAnnotation> linxSvAnnotations = List.of(svAnnotation());
        List<LinxBreakend> breakends = createBreakends();

        LinxBreakendInterpreter interpreter = createInterpreter(structuralVariants, linxSvAnnotations);
        com.hartwig.hmftools.datamodel.linx.LinxBreakend left = interpreter.interpret(breakends.get(0));
        assertNotNull(left);
        assertEquals("7", left.chromosome());
        assertEquals("q34", left.chromosomeBand());
        assertEquals(LinxBreakendType.BND, left.type());
        assertEquals(1, left.orientation());
        assertEquals(1.5D, left.junctionCopyNumber(), EPSILON);

        com.hartwig.hmftools.datamodel.linx.LinxBreakend right = interpreter.interpret(breakends.get(1));
        assertFallbackBreakendAttributes(right);
    }

    @NotNull
    private StructuralVariant createStructuralVariant(@NotNull StructuralVariantType type)
    {
        return PurpleTestUtils.createStructuralVariant("7", 50, "7", 200, type, 1.0, 1.0)
                .id("vcfId")
                .build();
    }

    @NotNull
    private List<LinxBreakend> createBreakends()
    {
        final LinxBreakend left = LinxTestFactory.breakendBuilder()
                .id(1)
                .reportedDisruption(true)
                .gene("gene")
                .isStart(true)
                .svId(1)
                .build();

        LinxBreakend right = LinxTestFactory.breakendBuilder()
                .id(2)
                .reportedDisruption(true)
                .gene("gene")
                .isStart(false)
                .svId(1)
                .build();

        return Lists.newArrayList(left, right);
    }

    @NotNull
    public static LinxSvAnnotation svAnnotation()
    {
        return LinxTestFactory.svAnnotationBuilder()
                .svId(1)
                .vcfId("vcfId")
                .junctionCopyNumberMin(1.0D)
                .junctionCopyNumberMax(2.0D)
                .build();
    }

    @NotNull
    private static LinxBreakendInterpreter createInterpreter(@NotNull final List<StructuralVariant> structuralVariants,
            @NotNull final List<LinxSvAnnotation> linxSvAnnotations)
    {
        return new LinxBreakendInterpreter(structuralVariants, linxSvAnnotations, TestEnsemblDataCacheFactory.loadTestCache());
    }

    public static void assertFallbackBreakendAttributes(@Nullable com.hartwig.hmftools.datamodel.linx.LinxBreakend breakend)
    {
        assertNotNull(breakend);
        assertEquals("", breakend.chromosome());
        assertEquals("", breakend.chromosomeBand());
        assertEquals(LinxBreakendType.BND, breakend.type());
        assertEquals(0, breakend.orientation());
        assertEquals(0.0D, breakend.junctionCopyNumber(), 0.0D);
    }
}