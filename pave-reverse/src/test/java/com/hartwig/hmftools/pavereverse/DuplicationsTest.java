package com.hartwig.hmftools.pavereverse;

import org.junit.Test;

public final class DuplicationsTest extends ReversePaveTestBase
{
    @Test
    public void duplicationAtExonBoundary()
    {
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("ARID1A", "D641dup");
        checkSingleChange(variant, "G", "GGAT", "chr1", 26760855);
    }

    @Test
    public void duplicationAtExonBoundaryNegativeStrand()
    {
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("BRAF", "V238dup");
        checkSingleChange(variant, "G", "GTAC", "chr7", 140_801_557);
    }

    @Test
    public void duplicationInRegionOfRepeatingAminoAcids()
    {
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("ARID1A", "P20_P21dup");
        checkChanges(variant,
            basesChange("A", "ACCCGCC", "chr1", 26_696_447),
            basesChange("C", "CCCGCCG", "chr1", 26_696_448),
            basesChange("G", "GCCGCCC", "chr1", 26_696_460)
        );
    }

    @Test
    public void duplicationOnReverseStrand()
    {
        BaseSequenceVariants v600 = reversePave.calculateProteinVariant("BRAF", "V600dup");
        checkSingleChange(v600, "T", "TCAC", "chr7", 140_753_334);

        BaseSequenceVariants t599 = reversePave.calculateProteinVariant("BRAF", "T599dup");
        // CAC TGT
        // CAC TGT TGT
        // CA CTG C TGT
        checkSingleChange(t599, "C", "CTGT", "chr7", 140_753_337);

        BaseSequenceVariants interval = reversePave.calculateProteinVariant("BRAF", "T599_V600dup");
        //     140_753_335
        //     |
        // K   V   T   A
        // TTT CAC TGT AGC
        // TTT CAC TGT AGC
        // TTT CAC TGT CAC TGT AGC
        checkSingleChange(interval, "T", "TTCACTG", "chr7", 140_753_333);
    }
}