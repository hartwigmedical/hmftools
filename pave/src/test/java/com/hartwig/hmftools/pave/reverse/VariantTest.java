package com.hartwig.hmftools.pave.reverse;

import java.util.List;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

public class VariantTest extends ReversePaveTestBase
{
    public static final String MAAQVAPAAS = "MAAQVAPAAS";
    private final int transcriptId = 111;
    private final ExonData exon0 = new ExonData(transcriptId, 9, 23, 1, 0, 0);
    private final ExonData exon1 = new ExonData(transcriptId, 27, 41, 2, 0, 0);
    GeneData gene = new GeneData("GeneId", "BLAH", "chr5", (byte) 1, 3, 36, "q13.1");
    TranscriptData transcript = new TranscriptData(500, "TrName", gene.GeneId, true, (byte) 1, 9, 42, 9, 42, "protein_coding");
    TranscriptAminoAcids taa = new TranscriptAminoAcids(gene.GeneId, gene.GeneName, transcript.TransName, true, MAAQVAPAAS);
    AminoAcidRange aar = new AminoAcidRange(aas(3, "A"), aas(4, "Q"));
    String exon0Bases = "ATGGCCGCGCAGGTC";
    String exon1Bases = "GCCCCCGCCGCCGCC";
    FixedStringGenome fsg = new FixedStringGenome("AAAGGGATC" + exon0Bases + "TTT" + exon1Bases + "TGATGATTT");

    public VariantTest()
    {
        transcript.setExons(List.of(exon0, exon1));
    }
}
