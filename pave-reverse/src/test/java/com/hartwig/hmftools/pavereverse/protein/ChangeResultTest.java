package com.hartwig.hmftools.pavereverse.protein;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.pavereverse.BaseSequenceChange;
import com.hartwig.hmftools.pavereverse.aminoacids.AminoAcidSequence;

import org.junit.Test;

public class ChangeResultTest extends VariantTest
{
    private final String Chr = "chr1";
    private final AminoAcidSequence mAminoAcidSequence = AminoAcidSequence.parse("MAAQVA");
    private final int Position = 13;
    private ChangeResult Change;
    private final String ExonBases = "ATGCCCCGCGCAGGTC";
    private final String AnchorBase = fsg.getBase(Chr, Position - 1);

    @Test
    public void toBaseSequenceSNVTest()
    {
        Change = new ChangeResult(mAminoAcidSequence, ExonBases, Position, "G", "A");
        checkChange(new BaseSequenceChange("G", "A", Chr, Position));
    }

    @Test
    public void toBaseSequenceDelTest()
    {
        Change = new ChangeResult(mAminoAcidSequence, ExonBases, Position, "GC", "");
        checkChange(new BaseSequenceChange(AnchorBase + "GC", AnchorBase, Chr, Position - 1));
    }

    @Test
    public void toBaseSequenceInsTest()
    {
        Change = new ChangeResult(mAminoAcidSequence, ExonBases, Position, "", "AT");
        checkChange(new BaseSequenceChange(AnchorBase, AnchorBase + "AT", Chr, Position - 1));
    }

    @Test
    public void toBaseSequenceDelInsTest()
    {
        Change = new ChangeResult(mAminoAcidSequence, ExonBases, Position, "GC", "ATA");
        checkChange(new BaseSequenceChange("GC", "ATA", Chr, Position));
    }

    private void checkChange(BaseSequenceChange expected)
    {
        assertEquals(expected, Change.asChange(Chr, fsg));
    }
}
