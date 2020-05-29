package com.hartwig.hmftools.serve.transvar;

import static org.junit.Assert.assertEquals;

import java.io.FileNotFoundException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;
import org.junit.Ignore;
import org.junit.Test;

public class TransvarInterpreterTest {

    private static final String REF_GENOME_FASTA_FILE = Resources.getResource("refgenome/ref.fasta").getPath();

    @Test
    public void canConvertSnvToHotspots() {
        TransvarRecord record = baseRecord().gdnaPosition(10)
                .gdnaRef("A")
                .gdnaAlt("C")
                .referenceCodon("TTA")
                .addCandidateCodons("GTA", "GTC", "GTG", "GTT")
                .build();

        List<VariantHotspot> hotspots = testInterpreter().convertRecordToHotspots(record, Strand.REVERSE);

        assertEquals(4, hotspots.size());

        assertHotspot(baseHotspot().position(10).ref("A").alt("C").build(), hotspots.get(0));
        assertHotspot(baseHotspot().position(8).ref("TAA").alt("GAC").build(), hotspots.get(1));
        assertHotspot(baseHotspot().position(8).ref("TAA").alt("CAC").build(), hotspots.get(2));
        assertHotspot(baseHotspot().position(8).ref("TAA").alt("AAC").build(), hotspots.get(3));
    }

    @Test
    public void canConvertMnvForwardStrandToHotspots() {
        TransvarRecord record = baseRecord().gdnaPosition(10)
                .gdnaRef("TA")
                .gdnaAlt("GC")
                .referenceCodon("TAC")
                .addCandidateCodons("GCA", "GCC", "GCG", "GCT")
                .build();

        List<VariantHotspot> hotspots = testInterpreter().convertRecordToHotspots(record, Strand.FORWARD);

        assertEquals(4, hotspots.size());

        assertHotspot(baseHotspot().position(10).ref("TAC").alt("GCA").build(), hotspots.get(0));
        assertHotspot(baseHotspot().position(10).ref("TA").alt("GC").build(), hotspots.get(1));
        assertHotspot(baseHotspot().position(10).ref("TAC").alt("GCG").build(), hotspots.get(2));
        assertHotspot(baseHotspot().position(10).ref("TAC").alt("GCT").build(), hotspots.get(3));
    }

    @Test
    public void canConvertMnvReverseStrandToHotspots() {
        TransvarRecord record = baseRecord().gdnaPosition(10)
                .gdnaRef("CA")
                .gdnaAlt("GC")
                .referenceCodon("TGG")
                .addCandidateCodons("GCA", "GCC", "GCG", "GCT")
                .build();

        List<VariantHotspot> hotspots = testInterpreter().convertRecordToHotspots(record, Strand.REVERSE);

        assertEquals(4, hotspots.size());

        assertHotspot(baseHotspot().position(9).ref("CCA").alt("TGC").build(), hotspots.get(0));
        assertHotspot(baseHotspot().position(9).ref("CCA").alt("GGC").build(), hotspots.get(1));
        assertHotspot(baseHotspot().position(10).ref("CA").alt("GC").build(), hotspots.get(2));
        assertHotspot(baseHotspot().position(9).ref("CCA").alt("AGC").build(), hotspots.get(3));
    }

    @Test
    public void canConvertDeletionToHotspots() {
        TransvarRecord record = baseRecord().gdnaPosition(5).gdnaRef("GAT").gdnaAlt("").indelLength(3).build();

        List<VariantHotspot> hotspots = testInterpreter().convertRecordToHotspots(record, Strand.FORWARD);

        assertEquals(1, hotspots.size());

        assertHotspot(baseHotspot().position(4).ref("CGAT").alt("C").build(), hotspots.get(0));
    }

    @Test
    public void canConvertInsertionToHotspots() {
        TransvarRecord forwardRecord = baseRecord().gdnaPosition(5).gdnaRef("").gdnaAlt("GAA").indelLength(3).build();

        List<VariantHotspot> forwardHotspots = testInterpreter().convertRecordToHotspots(forwardRecord, Strand.FORWARD);

        assertEquals(2, forwardHotspots.size());

        assertHotspot(baseHotspot().position(4).ref("C").alt("CGAG").build(), forwardHotspots.get(0));
        assertHotspot(baseHotspot().position(4).ref("C").alt("CGAA").build(), forwardHotspots.get(1));

        TransvarRecord reverseRecord = baseRecord().gdnaPosition(5).gdnaRef("").gdnaAlt("GAA").indelLength(3).build();

        List<VariantHotspot> reverseHotspots = testInterpreter().convertRecordToHotspots(reverseRecord, Strand.REVERSE);

        assertEquals(2, reverseHotspots.size());

        assertHotspot(baseHotspot().position(4).ref("C").alt("CAAA").build(), reverseHotspots.get(0));
        assertHotspot(baseHotspot().position(4).ref("C").alt("CGAA").build(), reverseHotspots.get(1));
    }

    @Test
    @Ignore
    public void canConvertCombinedDeletionsInsertionsToHotspot() {
        // TODO
        // See for example EGFR:p.L747_A750delinsP
    }

    @Test
    public void canConvertDuplicationToHotspot() {
        TransvarRecord record = baseRecord().gdnaPosition(5).gdnaRef("").gdnaAlt("").indelLength(3).build();

        List<VariantHotspot> hotspots = testInterpreter().convertRecordToHotspots(record, Strand.FORWARD);

        assertEquals(1, hotspots.size());

        assertHotspot(baseHotspot().position(4).ref("C").alt("CGAT").build(), hotspots.get(0));
    }

    private static void assertHotspot(@NotNull VariantHotspot expectedHotspot, @NotNull VariantHotspot actualHotspot) {
        assertEquals(expectedHotspot.chromosome(), actualHotspot.chromosome());
        assertEquals(expectedHotspot.position(), actualHotspot.position());
        assertEquals(expectedHotspot.ref(), actualHotspot.ref());
        assertEquals(expectedHotspot.alt(), actualHotspot.alt());
    }

    @NotNull
    private static ImmutableTransvarRecord.Builder baseRecord() {
        return ImmutableTransvarRecord.builder().transcript("irrelevant").chromosome("1");
    }

    @NotNull
    private static ImmutableVariantHotspotImpl.Builder baseHotspot() {
        return ImmutableVariantHotspotImpl.builder().chromosome("1");
    }

    @NotNull
    private static TransvarInterpreter testInterpreter() {
        try {
            return TransvarInterpreter.fromRefGenomeFastaFile(REF_GENOME_FASTA_FILE);
        } catch (FileNotFoundException exception) {
            throw new IllegalStateException("Cannot create test interpreter! Message=" + exception.getMessage());
        }
    }
}