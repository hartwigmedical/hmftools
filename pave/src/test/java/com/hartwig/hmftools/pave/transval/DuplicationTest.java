package com.hartwig.hmftools.pave.transval;

import java.util.List;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.junit.Before;
import org.junit.Test;

public class DuplicationTest extends TransvalTestBase
{
    private final int transcriptId = 111;
    private final ExonData exon0 = new ExonData(transcriptId, 9, 24, 1, 0, 0);
    private final ExonData exon1 = new ExonData(transcriptId, 27, 42, 2, 0, 0);
    GeneData gene = new GeneData("GeneId", "BLAH", "chr5", (byte) 1, 3, 36, "q13.1");
    TranscriptData transcript = new TranscriptData(500, "TrName", gene.GeneId, true, (byte) 1, 9, 42, 9, 42, "protein_coding");
    TranscriptAminoAcids taa = new TranscriptAminoAcids(gene.GeneId, gene.GeneName, transcript.TransName, true, "MAAQVAPAAS");
    AminoAcidRange aar = new AminoAcidRange(aas(3, "A"), aas(4, "Q"));
    FixedStringGenome fsg = new FixedStringGenome("AAAGGGATC"+ "ATGGCCGCGCAGGTC"  + "TTT" + "GCCCCCGCCGCCGCC" + "TGATGATTT");

    @Before
    public void setup()
    {
        transcript.setExons(List.of(exon0, exon1));
    }

    @Test
    public void calculateVariantTest()
    {
        Duplication duplication = new Duplication(gene, transcript, taa, aar);
        TransvalVariant variant = duplication.calculateVariant(fsg);
        checkSingleHotspot(variant, "C", "CGCGCAG", "chr5", 14);
    }

    private Duplication dup(String gene, String variant)
    {
        return transval.variationParser().parseDuplication(gene, variant);
    }
}
