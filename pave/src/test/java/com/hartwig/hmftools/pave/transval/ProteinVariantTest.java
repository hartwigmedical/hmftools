package com.hartwig.hmftools.pave.transval;

import java.util.Iterator;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.junit.Assert;
import org.junit.Test;

public class ProteinVariantTest extends VariantTest
{
    @Test
    public void calculateVariantTest()
    {
        final int iId = 111;
        GeneData g = new GeneData("GeneId", "BLAH", "chr5", (byte) 1, 10, 115, "q13.1");
        TranscriptData transcriptData = new TranscriptData(500, "TrName", g.GeneId, true, (byte) 1, 10, 115, 10, 115, "protein_coding");
        String aminoAcids = "LARGEFLAPPYELEPHANTSEARS";
        TranscriptAminoAcids transcriptAminoAcids = new TranscriptAminoAcids(g.GeneId, g.GeneName, transcriptData.TransName, true, aminoAcids);

        // TTAGCAAGGGGAGAATTCTTAGCACCACCATATGAATTAGAACCACACGCAAACACCTCAGAAGCAAGGTCA
        // L   A   R   G    E   F   L    A   P    P   Y   E    L   E   P    H   A   N    T   S   E    A   R   S
        // TTA|GCA|AGG|G GA|GAA|TTC|TT A|GCA|CCA| CCA|TAT|GAA| TTA|GAA|C CA|CAC|GCA|A AC|ACC|TCA|GA A|GCA|AGG|TCA
        // TTAGCAAGGG GAGAATTCTT AGCACCA CCATATGAA TTAGAAC CACACGCAA ACACCTCAGA AGCAAGGTCA
        // 10         10         7       9         7       9         10         10
        final String b0 = "TTAGCAAGGG";
        final ExonData e0 = new ExonData(iId, 10, 19, 1, 0, 1);
        final String b1 = "GAGAATTCTT";
        final ExonData e1 = new ExonData(iId, 22, 31, 2, 2, 2);
        final String b2 = "AGCACCA";
        final ExonData e2 = new ExonData(iId, 36, 42, 3, 1, 0);
        final String b3 = "CCATATGAA";
        final ExonData e3 = new ExonData(iId, 47, 55, 4, 0, 0);
        final String b4 = "TTAGAAC";
        final ExonData e4 = new ExonData(iId, 60, 66, 5, 0, 1);
        final String b5 = "CACACGCAA";
        final ExonData e5 = new ExonData(iId, 70, 78, 6, 2, 1);
        final String b6 = "ACACCTCAGA";
        final ExonData e6 = new ExonData(iId, 83, 92, 7, 2, 2);
        final String b7 = "AGCAAGGTCA";
        final ExonData e7 = new ExonData(iId, 97, 106, 8, 1, 0);
        // AAAGAGGATC TTAGCAAGGG TT GAGAATTCTT TTTT
        // AAAGAGGATCTTAGCAAGGG
        String allBases = "AAAGAGGATC" +
                b0 + "TT" +
                b1 + "AAAA" +
                b2 + "CCCC" +
                b3 + "GGGG" +
                b4 + "TTT" +
                b5 + "AAAA" +
                b6 + "CCCC" +
                b7 + "TGATGATTT";
        FixedStringGenome fixedGenome = new FixedStringGenome(allBases);
        transcriptData.setExons(List.of(e0, e1, e2, e3, e4, e5, e6, e7));

        ProteinVariant variant = new ProteinVariant(g, transcriptData, transcriptAminoAcids, 1, 2);
        //LARg gEFl lAP PYE LEp pHAn nTSe eARS
        checkExonReplacement("QQEFLAPPYELEPHANTSEARS", variant, 0, "QQ");
        checkExonReplacement("LARQQQAPPYELEPHANTSEARS", variant, 1, "QQQ");
        checkExonReplacement("LARGEFWWWPYELEPHANTSEARS", variant, 2, "WWW");
        checkExonReplacement("LARGEFLAPTTTLEPHANTSEARS", variant, 3, "TTT");
        checkExonReplacement("LARGEFLAPPYEWWWHANTSEARS", variant, 4, "WWW");
        checkExonReplacement("LARGEFLAPPYELERRRRTSEARS", variant, 5, "RRRR");
        checkExonReplacement("LARGEFLAPPYELEPHAFFFARS", variant, 6, "FFF");
        checkExonReplacement("LARGEFLAPPYELEPHANTSWWW", variant, 7, "WWW");
    }

    private void checkExonReplacement(String expected, ProteinVariant variant, int exon, String newAminoAcids)
    {
        AminoAcidSequence replacement = AminoAcidSequence.parse(newAminoAcids);
        AminoAcidSequence expectedAAs = AminoAcidSequence.parse(expected);
        Assert.assertEquals(expectedAAs, variant.replaceExonAminoAcids(exon, replacement));
    }
}
