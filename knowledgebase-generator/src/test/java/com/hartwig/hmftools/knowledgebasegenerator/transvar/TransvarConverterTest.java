package com.hartwig.hmftools.knowledgebasegenerator.transvar;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class TransvarConverterTest {

    @Test
    public void canConvertTransvarLineToRecord() {
        String line =
                "MTOR:p.L2230V\tENST00000361445 (protein_coding)\tMTOR\t-\tchr1:g.11182158A>C/c.6688T>G/p.L2230V\tinside_[cds_in_exon_48]"
                        + "\tCSQN=Missense;reference_codon=TTA;candidate_codons=GTA,GTC,GTG,GTT;candidate_mnv_variants="
                        + "chr1:g.11182156_11182158delTAAinsGAC,chr1:g.11182156_11182158delTAAinsCAC,chr1:g.11182156_11182158delTAAinsAAC;"
                        + "aliases=ENSP00000354558;source=Ensembl";

        TransvarRecord record = TransvarConverter.toTransvarRecord(line);

        assertEquals("ENST00000361445", record.transcript());
        assertEquals("1", record.chromosome());
        assertEquals(11182158, record.gdnaPosition());
        assertEquals("A", record.gdnaRef());
        assertEquals("C", record.gdnaAlt());
        assertEquals("TTA", record.referenceCodon());
        assertEquals("GTA", record.candidateCodons().get(0));
        assertEquals("GTC", record.candidateCodons().get(1));
        assertEquals("GTG", record.candidateCodons().get(2));
        assertEquals("GTT", record.candidateCodons().get(3));
    }

    @Test
    public void canConvertRecordToHotspots() {
        TransvarRecord record = ImmutableTransvarRecord.builder()
                .transcript("Irrelevant")
                .chromosome("1")
                .gdnaPosition(11182158)
                .gdnaRef("A")
                .gdnaAlt("C")
                .referenceCodon("TTA")
                .addCandidateCodons("GTA", "GTC", "GTG", "GTT")
                .build();

        List<VariantHotspot> hotspots = TransvarConverter.convertRecordToHotspots(record, Strand.REVERSE);

        assertEquals(4, hotspots.size());

        assertHotspot(chr1().position(11182158).ref("A").alt("C").build(), hotspots.get(0));
        assertHotspot(chr1().position(11182156).ref("TAA").alt("GAC").build(), hotspots.get(1));
        assertHotspot(chr1().position(11182156).ref("TAA").alt("CAC").build(), hotspots.get(2));
        assertHotspot(chr1().position(11182156).ref("TAA").alt("AAC").build(), hotspots.get(3));
    }

    private static void assertHotspot(@NotNull VariantHotspot expectedHotspot, @NotNull VariantHotspot actualHotspot) {
        assertEquals(expectedHotspot.chromosome(), actualHotspot.chromosome());
        assertEquals(expectedHotspot.position(), actualHotspot.position());
        assertEquals(expectedHotspot.ref(), actualHotspot.ref());
        assertEquals(expectedHotspot.alt(), actualHotspot.alt());
    }

    @NotNull
    private static ImmutableVariantHotspotImpl.Builder chr1() {
        return ImmutableVariantHotspotImpl.builder().chromosome("1");
    }
}