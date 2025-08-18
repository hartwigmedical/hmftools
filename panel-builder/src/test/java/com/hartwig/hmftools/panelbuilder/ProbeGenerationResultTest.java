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
        ProbeTarget target = ProbeTarget.exactRegion(region);
        String sequence = MockRefGenome.generateRandomBases(region.baseLength());
        return new Probe(target, sequence, metadata, null, null, 0, 0)
                .withEvalCriteria(new ProbeEvaluator.Criteria(1.0, 0.5, 0.1))
                .withRejectionReason(null);
    }

    @Test
    public void testAdd()
    {
        ProbeGenerationResult result1 = new ProbeGenerationResult(
                List.of(probe(new ChrBaseRegion("1", 10, 20), new TargetMetadata(TargetMetadata.Type.GENE, "1"))),
                List.of(new TargetRegion(new ChrBaseRegion("2", 20, 30), new TargetMetadata(TargetMetadata.Type.CUSTOM, "2"))),
                List.of(new TargetRegion(new ChrBaseRegion("3", 30, 40), new TargetMetadata(TargetMetadata.Type.CN_BACKBONE, "3"))),
                List.of(new RejectedRegion(new ChrBaseRegion("4", 40, 50), new TargetMetadata(TargetMetadata.Type.CUSTOM, "4"), "4"))
        );
        ProbeGenerationResult result2 = new ProbeGenerationResult(
                List.of(probe(new ChrBaseRegion("5", 50, 60), new TargetMetadata(TargetMetadata.Type.CN_BACKBONE, "5"))),
                List.of(new TargetRegion(new ChrBaseRegion("6", 60, 70), new TargetMetadata(TargetMetadata.Type.GENE, "6"))),
                List.of(new TargetRegion(new ChrBaseRegion("7", 70, 80), new TargetMetadata(TargetMetadata.Type.CUSTOM, "7"))),
                List.of(new RejectedRegion(new ChrBaseRegion("8", 80, 90), new TargetMetadata(TargetMetadata.Type.CN_BACKBONE, "8"), "8"))
        );
        ProbeGenerationResult expected = new ProbeGenerationResult(
                List.of(
                        probe(new ChrBaseRegion("1", 10, 20), new TargetMetadata(TargetMetadata.Type.GENE, "1")),
                        probe(new ChrBaseRegion("5", 50, 60), new TargetMetadata(TargetMetadata.Type.CN_BACKBONE, "5"))
                ),
                List.of(
                        new TargetRegion(new ChrBaseRegion("2", 20, 30), new TargetMetadata(TargetMetadata.Type.CUSTOM, "2")),
                        new TargetRegion(new ChrBaseRegion("6", 60, 70), new TargetMetadata(TargetMetadata.Type.GENE, "6"))
                ),
                List.of(
                        new TargetRegion(new ChrBaseRegion("3", 30, 40), new TargetMetadata(TargetMetadata.Type.CN_BACKBONE, "3")),
                        new TargetRegion(new ChrBaseRegion("7", 70, 80), new TargetMetadata(TargetMetadata.Type.CUSTOM, "7"))
                ),
                List.of(
                        new RejectedRegion(new ChrBaseRegion("4", 40, 50), new TargetMetadata(TargetMetadata.Type.CUSTOM, "4"), "4"),
                        new RejectedRegion(new ChrBaseRegion("8", 80, 90), new TargetMetadata(TargetMetadata.Type.CN_BACKBONE, "8"), "8")
                )
        );
        ProbeGenerationResult actual = result1.add(result2);
        assertEquals(expected, actual);
    }

    // TODO: test others

    @Test
    public void testRejectTarget()
    {
        TargetMetadata metadata = new TargetMetadata(TargetMetadata.Type.CUSTOM, "extra");
        TargetRegion target = new TargetRegion(
                new ChrBaseRegion("1", 10, 20),
                metadata);
        String reason = "rejected";
        ProbeGenerationResult actual = ProbeGenerationResult.rejectTarget(target, reason);
        ProbeGenerationResult expected = new ProbeGenerationResult(
                emptyList(),
                List.of(target),
                emptyList(),
                List.of(new RejectedRegion(target.region(), metadata, reason))
        );
        assertEquals(expected, actual);
    }
}
