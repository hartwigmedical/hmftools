package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

// Verifies the RNA output writers handle a multi-region (spliced) probe, which the DNA probe format cannot represent.
public class OutputWriterRnaTest
{
    private static final ProbeEvaluator.Criteria CRITERIA = new ProbeEvaluator.Criteria(0.5, 0.5, 1.0);

    @Test
    public void testWriteRnaSplicedProbe() throws IOException
    {
        // A spliced probe across two exons, each half the probe length, so the total is exactly the probe length.
        int half = PROBE_LENGTH / 2;
        ChrBaseRegion region1 = new ChrBaseRegion("1", 1001, 1000 + half);
        ChrBaseRegion region2 = new ChrBaseRegion("1", 2001, 2000 + half);
        SequenceDefinition definition = SequenceDefinition.spliced(List.of(region1, region2));
        Probe probe = new Probe(definition, new TargetedRange(0, PROBE_LENGTH), new TargetMetadata(TargetMetadata.Type.GENE_RNA, "GENE1"))
                .withSequence("A".repeat(PROBE_LENGTH))
                .withEvaluationCriteria(CRITERIA)
                .withEvaluationResult(EvaluationResult.accept())
                .withQualityScore(1.0)
                .withGcContent(0.5);

        Path outputDir = Files.createTempDirectory("panelbuilder-rna-output");
        try(OutputWriter writer = new OutputWriter(outputDir.toString(), null, false, true))
        {
            writer.rnaPanelOutput().writeProbes(List.of(probe));
        }

        List<String> tsvLines = Files.readAllLines(outputDir.resolve("rna_probes.tsv"));
        assertEquals(2, tsvLines.size());
        String dataRow = tsvLines.get(1);
        // Each segment is serialized as region:orientation (1/-1), joined by ';'. Both ref segments are genome-forward.
        assertTrue(dataRow.contains(region1 + ":1;" + region2 + ":1"));
        assertTrue(dataRow.contains("GENE_RNA"));

        // Each region of the spliced probe appears as its own BED row.
        List<String> bedLines = Files.readAllLines(outputDir.resolve("rna_probes.bed"));
        assertEquals(2, bedLines.size());
    }
}
