package com.hartwig.hmftools.pavereverse;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import com.hartwig.hmftools.pavereverse.aminoacids.AminoAcidSequence;

import org.junit.Test;

public final class InsertionsTest extends ReversePaveTestBase
{
    @Test
    public void insertionAtStartOfExon()
    {
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("VHL", "A5_E6insW");

        //RAE: AGG GCG GAG, 2nd G of the A is at 10_141_862
        checkSingleChange(variant, "G", "GTGG", "chr3", 10_141_862);

        // Same location but multiple codons for the inserted amino acid.
        variant = reversePave.calculateProteinVariant("VHL", "A5_E6insF");
        checkChanges(variant,
                basesChange("G", "GTTT", "chr3", 10_141_862),
                basesChange("G", "GTTC", "chr3", 10_141_862)
        );

        // Same location but inserting multiple amino acid.
        variant = reversePave.calculateProteinVariant("VHL", "A5_E6insFARM");
        Set<BaseSequenceChange> hotspots = variant.changes();
        assertEquals(1, hotspots.size());
        BaseSequenceChange hotspot = hotspots.iterator().next();
        assertEquals("G", hotspot.Ref);
        String codons = hotspot.Alt.substring(1);
        assertEquals("FARM", AminoAcidSequence.fromNucleotides(codons).sequence());
        assertEquals(10_141_862, hotspot.Position);
    }

    @Test
    public void insertionEGFRTest()
    {
        // Bug that was found by comparing with Transvar.
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("EGFR", "H773_V774insSH");
        // ...P H C V... CCC CAC GTG TGC... the second C of the H is at 55,181,328
        Set<BaseSequenceChange> hotspots = variant.changes();
        assertEquals(1, hotspots.size());
        BaseSequenceChange hotspot = hotspots.iterator().next();
        assertEquals("C", hotspot.Ref);
        String codons = hotspot.Alt.substring(1);
        assertEquals("SH", AminoAcidSequence.fromNucleotides(codons).sequence());
        assertEquals(55_181_328, hotspot.Position);
    }

    @Test
    public void transvarInsertionTest()
    {
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("EGFR", "P772_H773insY");
        checkChanges(variant,
                basesChange("C", "CTAC", "chr7", 55_181_325),
                basesChange("C", "CTAT", "chr7", 55_181_325)
        );
    }

    @Test
    public void insertLongSequence()
    {
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("VHL", "A5_E6insGLVQVTGSSDNEYFYVDFREYE");
        Set<BaseSequenceChange> hotspots = variant.changes();
        assertEquals(1, hotspots.size());
        BaseSequenceChange hotspot = hotspots.iterator().next();
        assertEquals("G", hotspot.Ref);
        String codons = hotspot.Alt.substring(1);
        assertEquals("GLVQVTGSSDNEYFYVDFREYE", AminoAcidSequence.fromNucleotides(codons).sequence());
        assertEquals(10_141_862, hotspot.Position);
    }

    @Test
    public void insertionOnReverseStrand()
    {
        var variant = reversePave.calculateProteinVariant("BRAF", "H510_V511insW");
        checkSingleChange(variant, "C", "CCCA", "chr7", 140_777_075);
    }
}