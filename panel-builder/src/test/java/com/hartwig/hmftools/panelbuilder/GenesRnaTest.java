package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.GenesRna.RnaRegionType.CODING;
import static com.hartwig.hmftools.panelbuilder.GenesRna.RnaRegionType.UTR_3;
import static com.hartwig.hmftools.panelbuilder.GenesRna.RnaRegionType.UTR_5;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertSame;

import java.util.List;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.panelbuilder.GenesRna.GeneOptions;
import com.hartwig.hmftools.panelbuilder.GenesRna.GeneTargets;
import com.hartwig.hmftools.panelbuilder.GenesRna.RnaTarget;

import org.junit.Test;

public class GenesRnaTest
{
    private static final byte FORWARD = 1;
    private static final byte REVERSE = -1;

    // Single transcript, four exons, coding span 2050..3099. Each exon is classified whole: any coding base -> coding.
    //   exonA [1000,1099] -> probe-space [0,100):   fully noncoding, below the coding span (5' fwd / 3' rev)
    //   exonB [2000,2099] -> probe-space [100,200):  part-coding (coding from 2050), so classified fully coding
    //   exonC [3000,3099] -> probe-space [200,300):  fully coding
    //   exonD [4000,4099] -> probe-space [300,400):  fully noncoding, above the coding span (3' fwd / 5' rev)
    private static List<TranscriptData> testTranscripts(byte strand)
    {
        TranscriptData transcript = new TranscriptData(
                1, "ENST1", "GENE1", true, strand, 1000, 4099, 2050, 3099, "protein_coding", null);
        transcript.setExons(List.of(
                new ExonData(1, 1000, 1099, 1, -1, -1),
                new ExonData(1, 2000, 2099, 2, -1, -1),
                new ExonData(1, 3000, 3099, 3, -1, -1),
                new ExonData(1, 4000, 4099, 4, -1, -1)));
        return List.of(transcript);
    }

    @Test
    public void testForwardStrandTargets()
    {
        GeneData gene = new GeneData("GENE1", "GENE1", "1", FORWARD, 1000, 4099, "");
        GeneTargets result = GenesRna.createTargets(gene, testTranscripts(FORWARD), new GeneOptions(true, true));

        assertEquals(400, result.mapping().length());
        // Part-coding exonB is covered whole as coding ([100,200) is the entire exon). Lower-genome UTR is 5', upper-genome UTR is 3'.
        assertEquals(
                List.of(
                        new RnaTarget(UTR_5, 0, 100),
                        new RnaTarget(CODING, 100, 200),
                        new RnaTarget(CODING, 200, 300),
                        new RnaTarget(UTR_3, 300, 400)),
                result.targets());
    }

    @Test
    public void testReverseStrandTargets()
    {
        // Same genome layout but reverse strand: the 5'/3' UTR assignment is flipped.
        GeneData gene = new GeneData("GENE1", "GENE1", "1", REVERSE, 1000, 4099, "");
        GeneTargets result = GenesRna.createTargets(gene, testTranscripts(REVERSE), new GeneOptions(true, true));

        assertEquals(
                List.of(
                        new RnaTarget(UTR_3, 0, 100),
                        new RnaTarget(CODING, 100, 200),
                        new RnaTarget(CODING, 200, 300),
                        new RnaTarget(UTR_5, 300, 400)),
                result.targets());
    }

    @Test
    public void testUtrsExcluded()
    {
        // Coding is always included; with both UTR flags off only the coding exons remain (the fully noncoding exons are dropped).
        GeneData gene = new GeneData("GENE1", "GENE1", "1", FORWARD, 1000, 4099, "");
        GeneTargets result = GenesRna.createTargets(gene, testTranscripts(FORWARD), new GeneOptions(false, false));

        assertEquals(
                List.of(
                        new RnaTarget(CODING, 100, 200),
                        new RnaTarget(CODING, 200, 300)),
                result.targets());
    }

    // Two exons; the terminal exon straddles the coding start with only a short coding tail (100b < a probe), so it is folded whole into a
    // single UTR target. Its noncoding bulk lies below the coding span, so it is the below-span UTR (3' on the reverse strand, e.g. a reverse
    // gene's 3' terminal exon like KRAS; 5' on the forward strand). The other exon [4000,4200] is fully coding.
    //   exonUtr    [1000,3099] -> probe-space [0,2100):     coding 3000..3099 (100b), noncoding below -> folded whole-exon UTR
    //   exonCoding [4000,4200] -> probe-space [2100,2301):  fully coding
    private static List<TranscriptData> straddlingBoundaryTranscript(byte strand)
    {
        TranscriptData transcript = new TranscriptData(
                1, "ENST1", "GENE1", true, strand, 1000, 4200, 3000, 4200, "protein_coding", null);
        transcript.setExons(List.of(
                new ExonData(1, 1000, 3099, 1, -1, -1),
                new ExonData(1, 4000, 4200, 2, -1, -1)));
        return List.of(transcript);
    }

    @Test
    public void testStraddlingBoundaryExonForward()
    {
        // Forward strand: the straddling exon's UTR is below the coding span, so 5'.
        GeneData gene = new GeneData("GENE1", "GENE1", "1", FORWARD, 1000, 4200, "");
        GeneTargets result = GenesRna.createTargets(gene, straddlingBoundaryTranscript(FORWARD), new GeneOptions(true, true));

        assertEquals(
                List.of(
                        new RnaTarget(UTR_5, 0, 2100),
                        new RnaTarget(CODING, 2100, 2301)),
                result.targets());
    }

    @Test
    public void testStraddlingBoundaryExonReverse()
    {
        // Reverse strand: the same below-span UTR is 3' (the mislabel-as-5' regression case, e.g. KRAS).
        GeneData gene = new GeneData("GENE1", "GENE1", "1", REVERSE, 1000, 4200, "");
        GeneTargets result = GenesRna.createTargets(gene, straddlingBoundaryTranscript(REVERSE), new GeneOptions(true, true));

        assertEquals(
                List.of(
                        new RnaTarget(UTR_3, 0, 2100),
                        new RnaTarget(CODING, 2100, 2301)),
                result.targets());
    }

    // Single exon [1000,1499] (500b), coding span 1200..1349: 200b of 5' UTR, 150b of coding, 150b of 3' UTR - each at least a probe length,
    // so the exon is partially coding and split into three separate targets.
    private static List<TranscriptData> partlyCodingTranscript(byte strand)
    {
        TranscriptData transcript = new TranscriptData(
                1, "ENST1", "GENE1", true, strand, 1000, 1499, 1200, 1349, "protein_coding", null);
        transcript.setExons(List.of(new ExonData(1, 1000, 1499, 1, -1, -1)));
        return List.of(transcript);
    }

    @Test
    public void testPartiallyCodingExonSplit()
    {
        GeneData gene = new GeneData("GENE1", "GENE1", "1", FORWARD, 1000, 1499, "");
        GeneTargets result = GenesRna.createTargets(gene, partlyCodingTranscript(FORWARD), new GeneOptions(true, true));

        // The coding part is its own target; the flanking UTR parts are separate (5' below the coding span, 3' above, on the forward strand).
        assertEquals(
                List.of(
                        new RnaTarget(UTR_5, 0, 200),
                        new RnaTarget(CODING, 200, 350),
                        new RnaTarget(UTR_3, 350, 500)),
                result.targets());
    }

    @Test
    public void testPartiallyCodingExonSplitUtrsExcluded()
    {
        // With UTRs disabled only the coding part of the split exon is covered - not the whole exon.
        GeneData gene = new GeneData("GENE1", "GENE1", "1", FORWARD, 1000, 1499, "");
        GeneTargets result = GenesRna.createTargets(gene, partlyCodingTranscript(FORWARD), new GeneOptions(false, false));

        assertEquals(List.of(new RnaTarget(CODING, 200, 350)), result.targets());
    }

    // Noncoding gene: a single noncoding transcript (no coding span) with two exons.
    private static List<TranscriptData> noncodingTranscript()
    {
        TranscriptData transcript = new TranscriptData(
                1, "ENST1", "GENE1", true, FORWARD, 1000, 1299, null, null, "lincRNA", null);
        transcript.setExons(List.of(
                new ExonData(1, 1000, 1099, 1, -1, -1),
                new ExonData(1, 1200, 1299, 2, -1, -1)));
        return List.of(transcript);
    }

    @Test
    public void testNoncodingGeneCoveredAsUtr()
    {
        // No coding span, so no 5'/3' distinction: each exon is one UTR target (labelled 5' when both UTRs are enabled).
        GeneData gene = new GeneData("GENE1", "GENE1", "1", FORWARD, 1000, 1299, "");
        GeneTargets result = GenesRna.createTargets(gene, noncodingTranscript(), new GeneOptions(true, true));

        assertEquals(
                List.of(new RnaTarget(UTR_5, 0, 100), new RnaTarget(UTR_5, 100, 200)),
                result.targets());
    }

    @Test
    public void testNoncodingGeneUtr3Only()
    {
        // With only 3' UTR enabled the noncoding exons are still covered (labelled 3').
        GeneData gene = new GeneData("GENE1", "GENE1", "1", FORWARD, 1000, 1299, "");
        GeneTargets result = GenesRna.createTargets(gene, noncodingTranscript(), new GeneOptions(false, true));

        assertEquals(
                List.of(new RnaTarget(UTR_3, 0, 100), new RnaTarget(UTR_3, 100, 200)),
                result.targets());
    }

    @Test
    public void testNoncodingGeneNoUtrRequested()
    {
        // Noncoding gene with no UTR requested: nothing to cover.
        GeneData gene = new GeneData("GENE1", "GENE1", "1", FORWARD, 1000, 1299, "");
        GeneTargets result = GenesRna.createTargets(gene, noncodingTranscript(), new GeneOptions(false, false));

        assertEquals(List.of(), result.targets());
    }

    @Test
    public void testParseTranscripts()
    {
        // Plain comma-separated.
        assertEquals(List.of("ENST1", "ENST2"), GenesRna.parseTranscripts("ENST1,ENST2"));
        // Excel wraps a comma-containing field in double quotes - the enclosing quotes must be stripped.
        assertEquals(List.of("ENST1", "ENST2"), GenesRna.parseTranscripts("\"ENST1,ENST2\""));
        // Single value, and whitespace around values.
        assertEquals(List.of("ENST1"), GenesRna.parseTranscripts("ENST1"));
        assertEquals(List.of("ENST1", "ENST2"), GenesRna.parseTranscripts("\"ENST1, ENST2\""));
        // Empty / null -> no transcripts.
        assertEquals(List.of(), GenesRna.parseTranscripts(null));
        assertEquals(List.of(), GenesRna.parseTranscripts(""));
    }

    @Test
    public void testResolveTranscript()
    {
        GeneData gene = new GeneData("GENE1", "GENE1", "1", FORWARD, 1000, 3099, "");
        TranscriptData enst1 = transcript("ENST1");
        TranscriptData enst2 = transcript("ENST2");
        List<TranscriptData> all = List.of(enst1, enst2);

        // Ensembl name match.
        assertSame(enst2, GenesRna.resolveTranscript(gene, all, "ENST2"));
        // Unknown Ensembl name.
        assertNull(GenesRna.resolveTranscript(gene, all, "ENST9"));
        // Only Ensembl ids are supported; a non-Ensembl (e.g. RefSeq) id is not matched.
        assertNull(GenesRna.resolveTranscript(gene, all, "NM_1"));
    }

    private static TranscriptData transcript(final String transName)
    {
        return new TranscriptData(1, transName, "GENE1", false, FORWARD, 1000, 2000, 1000, 2000, "protein_coding", null);
    }
}
