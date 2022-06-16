package com.hartwig.hmftools.common.linx;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.common.sv.linx.ImmutableLinxSvAnnotation;
import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.linx.LinxGermlineSv;
import com.hartwig.hmftools.common.sv.linx.LinxSvAnnotation;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GeneDisruptionFactoryTest {

    private static final double EPSILON = 1.0E-10;

    @Test
    public void canDetermineClusterId() {
        LinxBreakend breakend = LinxTestFactory.breakendBuilder()
                .svId(1)
                .gene("ROPN1B")
                .chromosome("3")
                .chrBand("p12")
                .type("INV")
                .junctionCopyNumber(1.12)
                .build();

        LinxSvAnnotation match = createTestAnnotationBuilder().svId(1).clusterId(2).build();
        assertEquals(2, (int) GeneDisruptionFactory.determineClusterId(Lists.newArrayList(match), breakend));

        LinxSvAnnotation noMatch = createTestAnnotationBuilder().svId(2).clusterId(1).build();
        assertNull(GeneDisruptionFactory.determineClusterId(Lists.newArrayList(noMatch), breakend));
    }

    @Test
    public void canConvertPairedBreakend() {
        ImmutableLinxBreakend.Builder pairedBreakendBuilder = LinxTestFactory.breakendBuilder()
                .svId(1)
                .gene("ROPN1B")
                .chromosome("3")
                .chrBand("p12")
                .type("INV")
                .junctionCopyNumber(1.12);
        List<LinxSvAnnotation> structuralVariants =
                Lists.newArrayList(createTestAnnotationBuilder().svId(1).clusterId(2).geneStart("ROPN1B").geneEnd("ROPN1B").build());

        List<LinxBreakend> pairedBreakends =
                Lists.newArrayList(pairedBreakendBuilder.exonUp(3).exonDown(4).undisruptedCopyNumber(4.3).build(),
                        pairedBreakendBuilder.exonUp(8).exonDown(9).undisruptedCopyNumber(2.1).build());

        List<GeneDisruption> geneDisruptions = GeneDisruptionFactory.convert(pairedBreakends, structuralVariants);

        assertEquals(1, geneDisruptions.size());

        GeneDisruption disruption = geneDisruptions.get(0);
        assertEquals("INV", disruption.type());
        assertEquals("3p12", disruption.location());
        assertEquals("ROPN1B", disruption.gene());
        assertEquals("Intron 3 -> Intron 8", disruption.range());
        assertEquals(3, disruption.firstAffectedExon());
        assertEquals(2.1, disruption.undisruptedCopyNumber(), EPSILON);

        Double copyNumber = disruption.junctionCopyNumber();
        assertNotNull(copyNumber);
        assertEquals(1.12, copyNumber, EPSILON);
    }

    @Test
    public void doesNotPairBreakendsOnDifferentGenes() {
        ImmutableLinxBreakend.Builder pairedBreakendBuilder = LinxTestFactory.breakendBuilder().svId(1);
        List<LinxSvAnnotation> structuralVariants =
                Lists.newArrayList(createTestAnnotationBuilder().svId(1).clusterId(2).geneStart("ROPN1B").geneEnd("ROPN1B").build());
        List<LinxBreakend> pairedDisruptions =
                Lists.newArrayList(pairedBreakendBuilder.gene("ROPN1B").svId(1).junctionCopyNumber(1.0).undisruptedCopyNumber(1.0).build(),
                        pairedBreakendBuilder.gene("SETD2").svId(1).junctionCopyNumber(1.0).undisruptedCopyNumber(2.3).build(),
                        pairedBreakendBuilder.gene("SETD2").svId(1).junctionCopyNumber(1.0).undisruptedCopyNumber(1.7).build());

        List<GeneDisruption> disruptions = GeneDisruptionFactory.convert(pairedDisruptions, structuralVariants);

        assertEquals(2, disruptions.size());
    }

    @NotNull
    private static ImmutableLinxSvAnnotation.Builder createTestAnnotationBuilder() {
        return ImmutableLinxSvAnnotation.builder()
                .vcfId(Strings.EMPTY)
                .svId(0)
                .clusterId(0)
                .clusterReason(Strings.EMPTY)
                .fragileSiteStart(true)
                .fragileSiteEnd(true)
                .isFoldback(true)
                .lineTypeStart(Strings.EMPTY)
                .lineTypeEnd(Strings.EMPTY)
                .junctionCopyNumberMin(0)
                .junctionCopyNumberMax(0)
                .geneStart(Strings.EMPTY)
                .geneEnd(Strings.EMPTY)
                .localTopologyIdStart(0)
                .localTopologyIdEnd(0)
                .localTopologyStart(Strings.EMPTY)
                .localTopologyEnd(Strings.EMPTY)
                .localTICountStart(0)
                .localTICountEnd(0);
    }
}