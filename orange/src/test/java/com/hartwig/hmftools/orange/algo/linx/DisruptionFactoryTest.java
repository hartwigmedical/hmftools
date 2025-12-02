package com.hartwig.hmftools.orange.algo.linx;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.hartwig.hmftools.datamodel.finding.Disruption;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.datamodel.linx.LinxGeneOrientation;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;

import org.junit.Test;

public class DisruptionFactoryTest
{
    @Test
    public void canDetermineClusterId()
    {
        LinxBreakend breakend = LinxOrangeTestFactory.breakendBuilder().svId(1).build();

        LinxSvAnnotation match = LinxOrangeTestFactory.svAnnotationBuilder().svId(1).clusterId(2).build();
        assertEquals(Integer.valueOf(2), DisruptionFactory.determineClusterId(List.of(match), breakend));

        LinxSvAnnotation noMatch = LinxOrangeTestFactory.svAnnotationBuilder().svId(2).clusterId(1).build();
        assertNull(DisruptionFactory.determineClusterId(List.of(noMatch), breakend));
    }

    @Test
    public void canConvertPairedBreakend()
    {
        ImmutableLinxBreakend.Builder pairedBreakendBuilder = LinxOrangeTestFactory.breakendBuilder()
                .svId(1)
                .gene("ROPN1B")
                .transcript("ENST1")
                .chromosome("3")
                .chromosomeBand("p12")
                .type(LinxBreakendType.INV)
                .junctionCopyNumber(1.12);
        List<LinxBreakend> pairedBreakends =
                List.of(pairedBreakendBuilder.exonUp(3).exonDown(4).undisruptedCopyNumber(4.3)
                                .geneOrientation(LinxGeneOrientation.UPSTREAM).build(),
                        pairedBreakendBuilder.exonUp(8).exonDown(9).undisruptedCopyNumber(2.1)
                                .geneOrientation(LinxGeneOrientation.DOWNSTREAM).build());

        List<Disruption> geneDisruptions = DisruptionFactory.createDisruptions(pairedBreakends, List.of(), true);

        assertEquals(1, geneDisruptions.size());

        Disruption disruption = geneDisruptions.get(0);
        assertEquals(LinxBreakendType.INV, disruption.type());
        //assertEquals("3p12", disruption.location());
        assertEquals("ROPN1B", disruption.gene());
        assertEquals(2.1, disruption.undisruptedCopies(), 0.000001);

        Double copyNumber = disruption.disruptedCopies();
        assertNotNull(copyNumber);
        assertEquals(1.12, copyNumber, 0.000001);
    }

    @Test
    public void upstreamBreakendOnly()
    {
        ImmutableLinxBreakend.Builder pairedBreakendBuilder = LinxOrangeTestFactory.breakendBuilder()
                .svId(1)
                .gene("ROPN1B")
                .transcript("ENST1")
                .chromosome("3")
                .chromosomeBand("p12")
                .type(LinxBreakendType.INV)
                .geneOrientation(LinxGeneOrientation.UPSTREAM)
                .junctionCopyNumber(1.12);
        List<LinxBreakend> pairedBreakends = List.of(pairedBreakendBuilder.exonUp(3).exonDown(4).undisruptedCopyNumber(4.3).build());

        List<Disruption> geneDisruptions = DisruptionFactory.createDisruptions(pairedBreakends, List.of(), true);

        assertEquals(1, geneDisruptions.size());

        Disruption disruption = geneDisruptions.get(0);
        assertEquals(LinxBreakendType.INV, disruption.type());
        //assertEquals("3p12", disruption.location());
        assertEquals("ROPN1B", disruption.gene());
        assertEquals(4.3, disruption.undisruptedCopies(), 0.000001);

        Double copyNumber = disruption.disruptedCopies();
        assertNotNull(copyNumber);
        assertEquals(1.12, copyNumber, 0.000001);
    }

    @Test
    public void doesNotPairBreakendsOnDifferentTranscripts()
    {
        ImmutableLinxBreakend.Builder pairedBreakendBuilder = LinxOrangeTestFactory.breakendBuilder().svId(1);
        List<LinxBreakend> pairedDisruptions = List.of(pairedBreakendBuilder.transcript("ENST 1").svId(1).build(),
                pairedBreakendBuilder.transcript("ENST 2").svId(1).build(),
                pairedBreakendBuilder.transcript("ENST 2").svId(1).build());

        List<Disruption> disruptions = DisruptionFactory.createDisruptions(pairedDisruptions, List.of(), true);

        assertEquals(2, disruptions.size());
    }
}
