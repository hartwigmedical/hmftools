package com.hartwig.hmftools.pavereverse;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public final class ReversePaveTest extends ReversePaveTestBase
{
    @Test
    public void calculateVariantsForSuppliedTranscriptId()
    {
        String transcriptId = "ENST00000361445";
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("MTOR", transcriptId, "L2230V");

        assertEquals(transcriptId, variant.transcriptName());
        assertEquals("1", variant.Chromosome);
        checkChanges(variant,
                basesChange("A", "C", "chr1", 11_122_101),
                basesChange("TAA", "GAC", "chr1", 11_122_099),
                basesChange("TAA", "CAC", "chr1", 11_122_099),
                basesChange("TAA", "AAC", "chr1", 11_122_099)
        );
    }

    @Test
    public void handleMultipleMatchingNonCanonicalTranscriptsThatReturnSameChanges()
    {
        BaseSequenceVariants variant = reversePave.calculateProteinVariantAllowingMultipleNonCanonicalTranscriptMatches("KIT", "N560_Y574del");
        checkChanges(variant,
                basesChange("TAAATGGAAACAATTATGTTTACATAGACCCAACACAACTTCCTTA", "T", "chr4", 54727456),
                basesChange("AAATGGAAACAATTATGTTTACATAGACCCAACACAACTTCCTTAT", "A", "chr4", 54727457)
        );
    }

    @Test
    public void handleVariantForWhichATranscriptIsIncomplete()
    {
        // One of the transcripts has a total length that is not a multiple of 3.
        BaseSequenceVariants variant = reversePave.calculateProteinVariantAllowingMultipleNonCanonicalTranscriptMatches("DYRK1A", "R437L");
        checkChanges(variant,
                basesChange("GA", "TT", "chr21", 37505353),
                basesChange("GA", "TG", "chr21", 37505353),
                basesChange("G", "T", "chr21", 37505353),
                basesChange("GA", "TC", "chr21", 37505353),
                basesChange("CG", "TT", "chr21", 37505352),
                basesChange("CGA", "TTG", "chr21", 37505352)
        );
    }

    @Test
    public void frameshiftForwardStrand()
    {
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("VHL", "Q132fs");
        // chr3, Q starts at 10146567
        // ...V N Q T...     ...GTT AAC CAA ACT...
        // AC>A at 565 gives ...GTT AA C AA ACT... which is ...V N K ....
        // Note that transvar gives CC>C @566
        checkSingleChange(variant, "AC", "A", "chr3", 10_146_565);

        // Similarly: VHL D9fs
        // ...N W D...     ...AAC TGG GAC G... the W starts at 10141869
        // TG>T @869 gives ...AAC TG GAC G...  which is ...N W T ....
        variant = reversePave.calculateProteinVariant("VHL", "D9fs");
        checkSingleChange(variant, "TG", "T", "chr3", 10_141_869);
    }

    @Test
    public void frameshiftReverseStrand()
    {
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("BRAF", "H585fs");
        checkChanges(variant,
                basesChange("CA", "C", "chr7", 140_753_379),
                basesChange("AT", "A", "chr7", 140_753_380),
                basesChange("TG", "T", "chr7", 140_753_381)
        );

        variant = reversePave.calculateProteinVariant("BRAF", "V600fs");
        checkChanges(variant,
                basesChange("CA", "C", "chr7", 140_753_335),
                basesChange("AC", "A", "chr7", 140_753_336)
        );
    }

    @Test
    public void stopGainedForwardStrand()
    {
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("VHL", "S65*");
        checkChanges(variant,
                basesChange("CG", "AA", "chr3", 10_142_041),
                basesChange("CG", "GA", "chr3", 10_142_041),
                basesChange("C", "A", "chr3", 10_142_041)
        );
    }

    @Test
    public void stopGainedReverseStrand()
    {
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("BRAF", "R603*");
        checkChanges(variant,
                basesChange("CG", "TA", "chr7", 140_753_327),
                basesChange("TCG", "CTA", "chr7", 140_753_326),
                basesChange("G", "A", "chr7", 140_753_328)
        );
    }

    @Test
    public void startLostForwardStrand()
    {
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("KIT", "M1?");
        checkChanges(variant,
                basesChange("A", "G", "chr4", 54_658_015),
                basesChange("A", "C", "chr4", 54_658_015),
                basesChange("A", "T", "chr4", 54_658_015),
                basesChange("T", "C", "chr4", 54_658_016),
                basesChange("T", "G", "chr4", 54_658_016),
                basesChange("T", "A", "chr4", 54_658_016),
                basesChange("G", "C", "chr4", 54_658_017),
                basesChange("G", "A", "chr4", 54_658_017),
                basesChange("G", "T", "chr4", 54_658_017)
        );
    }

    @Test
    public void startLostReverseStrand()
    {
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("BRCA1", "M1?");
        checkChanges(variant,
                basesChange("C", "T", "chr17", 43_124_094),
                basesChange("C", "G", "chr17", 43_124_094),
                basesChange("C", "A", "chr17", 43_124_094),
                basesChange("A", "C", "chr17", 43_124_095),
                basesChange("A", "G", "chr17", 43_124_095),
                basesChange("A", "T", "chr17", 43_124_095),
                basesChange("T", "A", "chr17", 43_124_096),
                basesChange("T", "C", "chr17", 43_124_096),
                basesChange("T", "G", "chr17", 43_124_096)
        );
    }
}