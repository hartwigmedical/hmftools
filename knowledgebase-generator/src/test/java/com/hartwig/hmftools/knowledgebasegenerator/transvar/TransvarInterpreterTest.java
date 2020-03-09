package com.hartwig.hmftools.knowledgebasegenerator.transvar;

import static org.junit.Assert.assertEquals;

import java.io.FileNotFoundException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class TransvarInterpreterTest {

    private static final String REF_GENOME_FASTA_FILE = Resources.getResource("refgenome/ref.fasta").getPath();

    @Test
    public void canConvertSNVRecordToHotspots() throws FileNotFoundException {
        TransvarRecord record = ImmutableTransvarRecord.builder()
                .transcript("Irrelevant")
                .chromosome("1")
                .gdnaPosition(10)
                .gdnaRef("A")
                .gdnaAlt("C")
                .referenceCodon("TTA")
                .addCandidateCodons("GTA", "GTC", "GTG", "GTT")
                .build();

        List<VariantHotspot> hotspots = testTransvarInterpreter().convertRecordToHotspots(record, Strand.REVERSE);

        assertEquals(4, hotspots.size());

        assertHotspot(chr1().position(10).ref("A").alt("C").build(), hotspots.get(0));
        assertHotspot(chr1().position(8).ref("TAA").alt("GAC").build(), hotspots.get(1));
        assertHotspot(chr1().position(8).ref("TAA").alt("CAC").build(), hotspots.get(2));
        assertHotspot(chr1().position(8).ref("TAA").alt("AAC").build(), hotspots.get(3));
    }

    @Test
    public void canConvertMNVRecordToHotspots() throws FileNotFoundException {
        TransvarRecord record = ImmutableTransvarRecord.builder()
                .transcript("Irrelevant")
                .chromosome("1")
                .gdnaPosition(10)
                .gdnaRef("TA")
                .gdnaAlt("GC")
                .referenceCodon("TAC")
                .addCandidateCodons("GCA", "GCC", "GCG", "GCT")
                .build();

        List<VariantHotspot> hotspots = testTransvarInterpreter().convertRecordToHotspots(record, Strand.FORWARD);

        assertEquals(4, hotspots.size());

        assertHotspot(chr1().position(10).ref("TAC").alt("GCA").build(), hotspots.get(0));
        assertHotspot(chr1().position(10).ref("TA").alt("GC").build(), hotspots.get(1));
        assertHotspot(chr1().position(10).ref("TAC").alt("GCG").build(), hotspots.get(2));
        assertHotspot(chr1().position(10).ref("TAC").alt("GCT").build(), hotspots.get(3));
    }

    @Test
    public void canConvertIndelRecordToHotspots() throws FileNotFoundException {
        TransvarRecord record = ImmutableTransvarRecord.builder()
                .transcript("Irrelevant")
                .chromosome("1")
                .gdnaPosition(5)
                .gdnaRef("GA")
                .gdnaAlt("")
                .build();

        List<VariantHotspot> hotspots = testTransvarInterpreter().convertRecordToHotspots(record, Strand.FORWARD);

        assertEquals(1, hotspots.size());

        assertHotspot(chr1().position(4).ref("CGA").alt("C").build(), hotspots.get(0));
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

    @NotNull
    private static TransvarInterpreter testTransvarInterpreter() throws FileNotFoundException {
        return TransvarInterpreter.fromRefGenomeFastaFile(REF_GENOME_FASTA_FILE);
    }

}