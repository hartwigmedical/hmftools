package com.hartwig.hmftools.pavereverse.protein;

import java.util.List;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.FixedStringGenome;
import com.hartwig.hmftools.pavereverse.aminoacids.AminoAcidRange;
import com.hartwig.hmftools.pavereverse.aminoacids.AminoAcidSequence;

import org.junit.Assert;
import org.junit.Test;

public class ProteinVariantTest extends VariantTest
{
    @Test
    public void replaceExonAminoAcids()
    {
        final int iId = 111;
        String aminoAcids = "MANGYFLAPPYELEPHANTSEARS";

        // M   A   N   G    Y   F   L    A   P    P   Y   E    L   E   P    H   A   N    T   S   E    A   R   S
        // ATG|GCA|AAC|G GA|TAT|TTC|TT A|GCA|CCA| CCA|TAT|GAA| TTA|GAA|C CA|CAC|GCA|A AC|ACC|TCA|GA A|GCA|AGG|TCA
        // ATGGCAAACG GATATTTCTT AGCACCA CCATATGAA TTAGAAC CACACGCAA ACACCTCAGA AGCAAGGTCA
        // 10         10         7       9         7       9         10         10
        final String b0 = "ATGGCAAACG";
        final ExonData e0 = new ExonData(iId, 10, 19, 1, 0, 1);
        final String b1 = "GATATTTCTT";
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
        int genomeLength = allBases.length();
        GeneData g = new GeneData("GeneId", "BLAH", "chr5", (byte) 1, 5, genomeLength - 5, "q13.1");
        TranscriptData transcriptData = new TranscriptData(500, "TrName", g.GeneId, true, (byte) 1, 10, genomeLength - 5, 10, genomeLength - 10, "protein_coding");
        TranscriptAminoAcids transcriptAminoAcids = new TranscriptAminoAcids(g.GeneId, g.GeneName, transcriptData.TransName, true, aminoAcids);

        FixedStringGenome fixedGenome = new FixedStringGenome(allBases);
        transcriptData.setExons(List.of(e0, e1, e2, e3, e4, e5, e6, e7));

        AminoAcidRange fsRange = new AminoAcidRange(aas(5, "F"), aas(5, "F"));
        ProteinVariant variant = new Frameshift(g, transcriptData, transcriptAminoAcids, fsRange);
        // If the replacement results in the loss of the initial M, return an empty sequence.
        AminoAcidSequence expectedAAs = variant.replaceExonAminoAcids(0, AminoAcidSequence.parse("QQ"));
//        Assert.assertEquals(0, expectedAAs.length());
        //MANg gYFl lAP PYE LEp pHAn nTSe eARS
        checkExonReplacement("MANQQQAPPYELEPHANTSEARS", variant, 1, "QQQ");
        checkExonReplacement("MANGYFWWWPYELEPHANTSEARS", variant, 2, "WWW");
        checkExonReplacement("MANGYFLAPTTTLEPHANTSEARS", variant, 3, "TTT");
        checkExonReplacement("MANGYFLAPPYEWWWHANTSEARS", variant, 4, "WWW");
        checkExonReplacement("MANGYFLAPPYELERRRRTSEARS", variant, 5, "RRRR");
        checkExonReplacement("MANGYFLAPPYELEPHAFFFARS", variant, 6, "FFF");
        checkExonReplacement("MANGYFLAPPYELEPHANTSWWW", variant, 7, "WWW");
    }

    private void checkExonReplacement(String expected, ProteinVariant variant, int exon, String newAminoAcids)
    {
        AminoAcidSequence replacement = AminoAcidSequence.parse(newAminoAcids);
        AminoAcidSequence expectedAAs = AminoAcidSequence.parse(expected);
        Assert.assertEquals(expectedAAs, variant.replaceExonAminoAcids(exon, replacement));
    }
}
