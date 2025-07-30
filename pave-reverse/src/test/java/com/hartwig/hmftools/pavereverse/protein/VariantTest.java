package com.hartwig.hmftools.pavereverse.protein;

import java.util.List;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.FixedStringGenome;
import com.hartwig.hmftools.pavereverse.ReversePaveTestBase;
import com.hartwig.hmftools.pavereverse.aminoacids.AminoAcidRange;

public class VariantTest extends ReversePaveTestBase
{
    public static final String MAAQVAPAAS = "MAAQVAPAAS";
    private final int transcriptId = 111;
    private final ExonData exon0 = new ExonData(transcriptId, 10, 24, 1, 0, 0);
    private final ExonData exon1 = new ExonData(transcriptId, 28, 42, 2, 0, 0);
    GeneData gene = new GeneData("GeneId", "BLAH", "chr5", (byte) 1, 3, 46, "q13.1");
    TranscriptData transcript = new TranscriptData(500, "TrName", gene.GeneId, true, (byte) 1, 9, 42, 10, 42, "protein_coding", "ref_seq_id");
    TranscriptAminoAcids taa = new TranscriptAminoAcids(gene.GeneId, gene.GeneName, transcript.TransName, true, MAAQVAPAAS);
    AminoAcidRange aar = new AminoAcidRange(aas(3, "A"), aas(4, "Q"));
    String exon0Bases = "ATGGCCGCGCAGGTC";
    String exon1Bases = "GCCCCCGCCGCCGCC";
    // 1        |10        20   |   |  30        40 |       50
    // 123456789|012345678901234|567|890123456789012|345678901
    // AAAGGGATC|ATGGCCGCGCAGGTC|TTT|GCCCCCGCCGCCGCC|TGATGATTT
    //          |M  A  A  Q  V  |   |A  P  A  A  S  |
    FixedStringGenome fsg = new FixedStringGenome("AAAGGGATC" + exon0Bases + "TTT" + exon1Bases + "TGATGATTT");

    public VariantTest()
    {
        transcript.setExons(List.of(exon0, exon1));
    }
}
