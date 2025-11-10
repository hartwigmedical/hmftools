package com.hartwig.hmftools.panelbuilder;

import static java.util.Collections.emptyList;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

public class ProbeGenerationResultTest
{
    public static Probe probe(final ChrBaseRegion region, final TargetMetadata metadata)
    {
        SequenceDefinition definition = SequenceDefinition.singleRegion(region);
        String sequence = MockRefGenome.generateRandomBases(region.baseLength());
        TargetedRange targetedRange = TargetedRange.wholeRegion(definition.baseLength());
        return new Probe(definition, sequence, targetedRange, metadata, null, null, 0.0, 0.0)
                .withEvalCriteria(new ProbeEvaluator.Criteria(1.0, 0.5, 0.1))
                .withRejectionReason(null);
    }

    @Test
    public void testAdd()
    {
        ProbeGenerationResult result1 = new ProbeGenerationResult(
                List.of(probe(new ChrBaseRegion("1", 10, 20), new TargetMetadata(TargetMetadata.Type.GENE, "1"))),
                List.of(new TargetRegion(new ChrBaseRegion("2", 20, 30), new TargetMetadata(TargetMetadata.Type.CUSTOM_REGION, "2"))),
                List.of(new RejectedRegion(new ChrBaseRegion("4", 40, 50), new TargetMetadata(TargetMetadata.Type.CUSTOM_REGION, "4")))
        );
        ProbeGenerationResult result2 = new ProbeGenerationResult(
                List.of(probe(new ChrBaseRegion("5", 50, 60), new TargetMetadata(TargetMetadata.Type.CN_BACKBONE, "5"))),
                List.of(new TargetRegion(new ChrBaseRegion("6", 60, 70), new TargetMetadata(TargetMetadata.Type.GENE, "6"))),
                List.of(new RejectedRegion(new ChrBaseRegion("8", 80, 90), new TargetMetadata(TargetMetadata.Type.CN_BACKBONE, "8")))
        );
        ProbeGenerationResult expected = new ProbeGenerationResult(
                List.of(
                        probe(new ChrBaseRegion("1", 10, 20), new TargetMetadata(TargetMetadata.Type.GENE, "1")),
                        probe(new ChrBaseRegion("5", 50, 60), new TargetMetadata(TargetMetadata.Type.CN_BACKBONE, "5"))
                ),
                List.of(
                        new TargetRegion(new ChrBaseRegion("2", 20, 30), new TargetMetadata(TargetMetadata.Type.CUSTOM_REGION, "2")),
                        new TargetRegion(new ChrBaseRegion("6", 60, 70), new TargetMetadata(TargetMetadata.Type.GENE, "6"))
                ),
                List.of(
                        new RejectedRegion(new ChrBaseRegion("4", 40, 50), new TargetMetadata(TargetMetadata.Type.CUSTOM_REGION, "4")),
                        new RejectedRegion(new ChrBaseRegion("8", 80, 90), new TargetMetadata(TargetMetadata.Type.CN_BACKBONE, "8"))
                )
        );
        ProbeGenerationResult actual = result1.add(result2);
        assertEquals(expected, actual);
    }

    @Test
    public void testAlreadyCoveredTarget()
    {
        TargetMetadata metadata = new TargetMetadata(TargetMetadata.Type.CUSTOM_REGION, "extra");
        TargetRegion target = new TargetRegion(
                new ChrBaseRegion("1", 10, 20),
                metadata);
        ProbeGenerationResult actual = ProbeGenerationResult.alreadyCoveredTargets(List.of(target));
        ProbeGenerationResult expected = new ProbeGenerationResult(
                emptyList(),
                List.of(target),
                emptyList()
        );
        assertEquals(expected, actual);
    }

    @Test
    public void testRejectTargets()
    {
        TargetMetadata metadata = new TargetMetadata(TargetMetadata.Type.CUSTOM_REGION, "extra");
        TargetRegion target = new TargetRegion(
                new ChrBaseRegion("1", 10, 20),
                metadata);
        ProbeGenerationResult actual = ProbeGenerationResult.rejectTargets(List.of(target));
        ProbeGenerationResult expected = new ProbeGenerationResult(
                emptyList(),
                List.of(target),
                List.of(new RejectedRegion(target.region(), metadata))
        );
        assertEquals(expected, actual);
    }
}
