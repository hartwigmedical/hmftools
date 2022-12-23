package com.hartwig.hmftools.orange.report.datamodel;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class BreakendEntryFactoryTest {

    private static final double EPSILON = 1.0E-10;

    @Test
    public void canCreateBreakendEntries() {
        LinxBreakend breakend = LinxTestFactory.breakendBuilder()
                .svId(1)
                .chromosome("1")
                .chrBand("p12.1")
                .gene("gene")
                .canonical(true)
                .exonUp(12)
                .exonDown(12)
                .geneOrientation("Upstream")
                .type(StructuralVariantType.DEL)
                .junctionCopyNumber(1.2)
                .undisruptedCopyNumber(1.4)
                .build();

        LinxSvAnnotation variant = LinxTestFactory.svAnnotationBuilder().svId(1).clusterId(2).build();

        List<BreakendEntry> entries = BreakendEntryFactory.create(Lists.newArrayList(breakend), Lists.newArrayList(variant));

        assertEquals(1, entries.size());

        BreakendEntry entry = entries.get(0);
        assertEquals("1p12.1", entry.location());
        assertEquals("gene", entry.gene());
        assertTrue(entry.canonical());
        assertEquals(12, entry.exonUp());
        assertEquals(StructuralVariantType.DEL, entry.type());
        assertEquals("Exon 12 Upstream", entry.range());
        assertEquals(2, entry.clusterId());
        assertEquals(1.2, entry.junctionCopyNumber(), EPSILON);
        assertEquals(1.4, entry.undisruptedCopyNumber(), EPSILON);
    }

    @Test (expected = IllegalStateException.class)
    public void crashOnMissingSvAnnotation() {
        LinxBreakend breakend = LinxTestFactory.breakendBuilder().svId(1).build();
        LinxSvAnnotation variant = LinxTestFactory.svAnnotationBuilder().svId(2).build();

        BreakendEntryFactory.create(Lists.newArrayList(breakend), Lists.newArrayList(variant));
    }

    @Test
    public void canGenerateRangeField() {
        assertEquals("Exon 4 Upstream", BreakendEntryFactory.range(create(4, 4, "Upstream")));
        assertEquals("Intron 4 Downstream", BreakendEntryFactory.range(create(4, 5, "Downstream")));
        assertEquals("Promoter Region Upstream", BreakendEntryFactory.range(create(0, 2, "Upstream")));
        assertEquals(Strings.EMPTY, BreakendEntryFactory.range(create(-1, -1, Strings.EMPTY)));
    }

    @NotNull
    private static LinxBreakend create(int exonUp, int exonDown, @NotNull String geneOrientation) {
        return LinxTestFactory.breakendBuilder().exonUp(exonUp).exonDown(exonDown).geneOrientation(geneOrientation).build();
    }
}