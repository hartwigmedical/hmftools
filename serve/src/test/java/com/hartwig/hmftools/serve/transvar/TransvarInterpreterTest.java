package com.hartwig.hmftools.serve.transvar;

import static com.hartwig.hmftools.serve.transvar.TransvarTestFactory.testInterpreter37;
import static com.hartwig.hmftools.serve.transvar.TransvarTestFactory.testInterpreter38;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarComplexInsertDelete;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarDeletion;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarDuplication;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarFrameshift;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarInsertion;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarRecord;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarSnvMnv;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarRecord;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class TransvarInterpreterTest {

    @Test
    public void canConvertSnvToHotspots() {
        TransvarRecord record = baseRecord().gdnaPosition(10)
                .annotation(ImmutableTransvarSnvMnv.builder()
                        .gdnaRef("A")
                        .gdnaAlt("C")
                        .referenceCodon("TTA")
                        .addCandidateCodons("GTA", "GTC", "GTG", "GTT")
                        .build())
                .build();

        List<VariantHotspot> hotspots = testInterpreter37().convertRecordToHotspots(record, Strand.REVERSE);

        assertEquals(4, hotspots.size());

        assertHotspot(baseHotspot().position(10).ref("A").alt("C").build(), hotspots.get(0));
        assertHotspot(baseHotspot().position(8).ref("TAA").alt("GAC").build(), hotspots.get(1));
        assertHotspot(baseHotspot().position(8).ref("TAA").alt("CAC").build(), hotspots.get(2));
        assertHotspot(baseHotspot().position(8).ref("TAA").alt("AAC").build(), hotspots.get(3));
    }

    @Test
    public void canConvertMnvForwardStrandToHotspots() {
        TransvarRecord record = baseRecord().gdnaPosition(10)
                .annotation(ImmutableTransvarSnvMnv.builder()
                        .gdnaRef("TA")
                        .gdnaAlt("GC")
                        .referenceCodon("TAC")
                        .addCandidateCodons("GCA", "GCC", "GCG", "GCT")
                        .build())
                .build();

        List<VariantHotspot> hotspots = testInterpreter37().convertRecordToHotspots(record, Strand.FORWARD);

        assertEquals(4, hotspots.size());

        assertHotspot(baseHotspot().position(10).ref("TAC").alt("GCA").build(), hotspots.get(0));
        assertHotspot(baseHotspot().position(10).ref("TA").alt("GC").build(), hotspots.get(1));
        assertHotspot(baseHotspot().position(10).ref("TAC").alt("GCG").build(), hotspots.get(2));
        assertHotspot(baseHotspot().position(10).ref("TAC").alt("GCT").build(), hotspots.get(3));
    }

    @Test
    public void canConvertMnvReverseStrandToHotspots() {
        TransvarRecord record = baseRecord().gdnaPosition(10)
                .annotation(ImmutableTransvarSnvMnv.builder()
                        .gdnaRef("CA")
                        .gdnaAlt("GC")
                        .referenceCodon("TGG")
                        .addCandidateCodons("GCA", "GCC", "GCG", "GCT")
                        .build())
                .build();

        List<VariantHotspot> hotspots = testInterpreter37().convertRecordToHotspots(record, Strand.REVERSE);

        assertEquals(4, hotspots.size());

        assertHotspot(baseHotspot().position(9).ref("CCA").alt("TGC").build(), hotspots.get(0));
        assertHotspot(baseHotspot().position(9).ref("CCA").alt("GGC").build(), hotspots.get(1));
        assertHotspot(baseHotspot().position(10).ref("CA").alt("GC").build(), hotspots.get(2));
        assertHotspot(baseHotspot().position(9).ref("CCA").alt("AGC").build(), hotspots.get(3));
    }

    @Test
    public void canInterpretSnvSpanningMultipleExons() {
        TransvarRecord record = baseRecord().variantSpanMultipleExons(true)
                .gdnaPosition(10)
                .annotation(ImmutableTransvarSnvMnv.builder()
                        .gdnaRef("A")
                        .gdnaAlt("C")
                        .referenceCodon("TTA")
                        .addCandidateCodons("GTA", "GTC", "GTG", "GTT")
                        .build())
                .build();

        List<VariantHotspot> hotspots = testInterpreter37().convertRecordToHotspots(record, Strand.REVERSE);

        assertEquals(1, hotspots.size());

        assertHotspot(baseHotspot().position(10).ref("A").alt("C").build(), hotspots.get(0));

        TransvarRecord record2 = baseRecord().variantSpanMultipleExons(true)
                .gdnaPosition(10)
                .annotation(ImmutableTransvarSnvMnv.builder()
                        .gdnaRef("TTA")
                        .gdnaAlt("GTA")
                        .referenceCodon("TTA")
                        .addCandidateCodons("GTA")
                        .build())
                .build();

        List<VariantHotspot> hotspots2 = testInterpreter37().convertRecordToHotspots(record2, Strand.FORWARD);

        assertEquals(0, hotspots2.size());
    }

    @Test
    public void canConvertDeletionToHotspots() {
        TransvarRecord record = baseRecord().gdnaPosition(5)
                .annotation(ImmutableTransvarDeletion.builder().deletedBaseCount(3).leftAlignedGDNAPosition(5).build())
                .build();

        List<VariantHotspot> hotspots = testInterpreter37().convertRecordToHotspots(record, Strand.FORWARD);

        assertEquals(1, hotspots.size());

        assertHotspot(baseHotspot().position(4).ref("CGAT").alt("C").build(), hotspots.get(0));
    }

    @Test
    public void canConvertDeletionWithMultipleAlignmentsToHotspots() {
        // In this situation the mutation is "GATCGATC -> GATC",
        //  Normally the unaligned DNA would be 1 here but that would imply we need to read the 0th ref base.
        TransvarRecord record = baseRecord().gdnaPosition(5)
                .annotation(ImmutableTransvarDeletion.builder().deletedBaseCount(4).leftAlignedGDNAPosition(2).build())
                .build();

        List<VariantHotspot> hotspots = testInterpreter37().convertRecordToHotspots(record, Strand.FORWARD);

        assertEquals(4, hotspots.size());

        assertHotspot(baseHotspot().position(1).ref("GATCG").alt("G").build(), hotspots.get(0));
        assertHotspot(baseHotspot().position(2).ref("ATCGA").alt("A").build(), hotspots.get(1));
        assertHotspot(baseHotspot().position(3).ref("TCGAT").alt("T").build(), hotspots.get(2));
        assertHotspot(baseHotspot().position(4).ref("CGATC").alt("C").build(), hotspots.get(3));
    }

    @Test
    public void canConvertInsertionsToHotspots() {
        TransvarRecord forwardRecord = baseRecord().gdnaPosition(5)
                .annotation(ImmutableTransvarInsertion.builder().insertedBases("GAA").leftAlignedGDNAPosition(4).build())
                .build();

        List<VariantHotspot> forwardHotspots = testInterpreter37().convertRecordToHotspots(forwardRecord, Strand.FORWARD);

        assertEquals(2, forwardHotspots.size());

        assertHotspot(baseHotspot().position(5).ref("G").alt("GGAG").build(), forwardHotspots.get(0));
        assertHotspot(baseHotspot().position(5).ref("G").alt("GGAA").build(), forwardHotspots.get(1));

        TransvarRecord reverseRecord = baseRecord().gdnaPosition(5)
                .annotation(ImmutableTransvarInsertion.builder().insertedBases("TGC").leftAlignedGDNAPosition(5).build())
                .build();

        List<VariantHotspot> reverseHotspots = testInterpreter37().convertRecordToHotspots(reverseRecord, Strand.REVERSE);

        assertEquals(4, reverseHotspots.size());

        assertHotspot(baseHotspot().position(5).ref("G").alt("GAGC").build(), reverseHotspots.get(0));
        assertHotspot(baseHotspot().position(5).ref("G").alt("GCGC").build(), reverseHotspots.get(1));
        assertHotspot(baseHotspot().position(5).ref("G").alt("GTGC").build(), reverseHotspots.get(2));
        assertHotspot(baseHotspot().position(5).ref("G").alt("GGGC").build(), reverseHotspots.get(3));
    }

    @Test
    public void canConvertComplexDeletionInsertionsToHotspots() {
        TransvarRecord oneAminoAcidInsert = baseRecord().gdnaPosition(2)
                .annotation(ImmutableTransvarComplexInsertDelete.builder()
                        .deletedBaseCount(6)
                        .insertedSequence("AAT")
                        .addCandidateAlternativeCodons("AAT")
                        .addCandidateAlternativeCodons("GAT")
                        .build())
                .build();

        List<VariantHotspot> hotspots1 = testInterpreter37().convertRecordToHotspots(oneAminoAcidInsert, Strand.FORWARD);

        assertEquals(2, hotspots1.size());

        assertHotspot(baseHotspot().position(2).ref("ATCG").alt("A").build(), hotspots1.get(0));
        assertHotspot(baseHotspot().position(2).ref("ATCG").alt("G").build(), hotspots1.get(1));

        TransvarRecord twoAminoAcidInsert = baseRecord().gdnaPosition(2)
                .annotation(ImmutableTransvarComplexInsertDelete.builder()
                        .deletedBaseCount(3)
                        .insertedSequence("AGGGTC")
                        .addCandidateAlternativeCodons("AGGGTC")
                        .addCandidateAlternativeCodons("ATTTTC")
                        .build())
                .build();

        List<VariantHotspot> hotspots2 = testInterpreter37().convertRecordToHotspots(twoAminoAcidInsert, Strand.FORWARD);

        assertEquals(1, hotspots2.size());

        assertHotspot(baseHotspot().position(2).ref("A").alt("AGGG").build(), hotspots2.get(0));
    }

    @Test
    public void canConvertComplexDeletionInsertionOnReverseStrand() {
        TransvarRecord oneAminoAcidInsert = baseRecord().gdnaPosition(2)
                .annotation(ImmutableTransvarComplexInsertDelete.builder()
                        .deletedBaseCount(3)
                        .insertedSequence("TAA")
                        .addCandidateAlternativeCodons("TTA")
                        .addCandidateAlternativeCodons("GCG")
                        .build())
                .build();

        List<VariantHotspot> hotspots1 = testInterpreter37().convertRecordToHotspots(oneAminoAcidInsert, Strand.REVERSE);

        assertEquals(2, hotspots1.size());

        assertHotspot(baseHotspot().position(2).ref("ATC").alt("TAA").build(), hotspots1.get(0));
        assertHotspot(baseHotspot().position(2).ref("AT").alt("CG").build(), hotspots1.get(1));
    }

    @Test
    public void canReduceComplexityOnComplexDelInsHotspots() {
        VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").position(10).ref("ATGTTA").alt("ATCCTA").build();

        VariantHotspot simplifiedHotspot = testInterpreter37().reduceComplexityForComplexInsDel(hotspot);

        assertEquals(12, simplifiedHotspot.position());
        assertEquals("GT", simplifiedHotspot.ref());
        assertEquals("CC", simplifiedHotspot.alt());
    }

    @Test
    public void canConvertDuplicationToHotspot37() {
        TransvarRecord record =
                baseRecord().gdnaPosition(5).annotation(ImmutableTransvarDuplication.builder().duplicatedBaseCount(3).build()).build();

        List<VariantHotspot> hotspots = testInterpreter37().convertRecordToHotspots(record, Strand.FORWARD);

        assertEquals(1, hotspots.size());

        assertHotspot(baseHotspot().position(4).ref("C").alt("CGAT").build(), hotspots.get(0));
    }

    @Test
    public void canConvertDuplicationToHotspot38() {
        TransvarRecord record =
                baseRecord().gdnaPosition(5).annotation(ImmutableTransvarDuplication.builder().duplicatedBaseCount(3).build()).build();

        List<VariantHotspot> hotspots = testInterpreter38().convertRecordToHotspots(record, Strand.FORWARD);

        assertEquals(1, hotspots.size());

        assertHotspot(ImmutableVariantHotspotImpl.builder().chromosome("chr1").position(4).ref("C").alt("CGAT").build(), hotspots.get(0));
    }

    @Test
    public void canConvertFrameshiftToHotspotOnForwardStrand() {
        TransvarRecord record = baseRecord().gdnaPosition(2)
                .annotation(ImmutableTransvarFrameshift.builder().isFrameshiftInsideStartCodon(false).build())
                .build();
        List<VariantHotspot> hotspots = testInterpreter37().convertRecordToHotspots(record, Strand.FORWARD);

        assertEquals(10, hotspots.size());
        assertHotspot(baseHotspot().position(2).ref("A").alt("AG").build(), hotspots.get(0));
        assertHotspot(baseHotspot().position(2).ref("A").alt("AA").build(), hotspots.get(1));
        assertHotspot(baseHotspot().position(2).ref("A").alt("AT").build(), hotspots.get(2));
        assertHotspot(baseHotspot().position(2).ref("A").alt("AC").build(), hotspots.get(3));

        assertHotspot(baseHotspot().position(3).ref("T").alt("TG").build(), hotspots.get(4));
        assertHotspot(baseHotspot().position(3).ref("T").alt("TA").build(), hotspots.get(5));
        assertHotspot(baseHotspot().position(3).ref("T").alt("TT").build(), hotspots.get(6));

        assertHotspot(baseHotspot().position(2).ref("AT").alt("A").build(), hotspots.get(7));

        assertHotspot(baseHotspot().position(2).ref("ATC").alt("A").build(), hotspots.get(8));
        assertHotspot(baseHotspot().position(3).ref("TCG").alt("T").build(), hotspots.get(9));
    }

    @Test
    public void canConvertFrameshiftToHotspotOnReverseStrand() {
        TransvarRecord record = baseRecord().gdnaPosition(6)
                .annotation(ImmutableTransvarFrameshift.builder().isFrameshiftInsideStartCodon(false).build())
                .build();
        List<VariantHotspot> hotspots = testInterpreter37().convertRecordToHotspots(record, Strand.REVERSE);

        assertEquals(10, hotspots.size());
        assertHotspot(baseHotspot().position(5).ref("G").alt("GA").build(), hotspots.get(0));
        assertHotspot(baseHotspot().position(5).ref("G").alt("GT").build(), hotspots.get(1));
        assertHotspot(baseHotspot().position(5).ref("G").alt("GC").build(), hotspots.get(2));

        assertHotspot(baseHotspot().position(6).ref("A").alt("AG").build(), hotspots.get(3));
        assertHotspot(baseHotspot().position(6).ref("A").alt("AA").build(), hotspots.get(4));
        assertHotspot(baseHotspot().position(6).ref("A").alt("AT").build(), hotspots.get(5));
        assertHotspot(baseHotspot().position(6).ref("A").alt("AC").build(), hotspots.get(6));

        assertHotspot(baseHotspot().position(5).ref("GA").alt("G").build(), hotspots.get(7));

        assertHotspot(baseHotspot().position(3).ref("TCG").alt("T").build(), hotspots.get(8));
        assertHotspot(baseHotspot().position(4).ref("CGA").alt("C").build(), hotspots.get(9));
    }

    private static void assertHotspot(@NotNull VariantHotspot expectedHotspot, @NotNull VariantHotspot actualHotspot) {
        assertEquals(expectedHotspot.chromosome(), actualHotspot.chromosome());
        assertEquals(expectedHotspot.position(), actualHotspot.position());
        assertEquals(expectedHotspot.ref(), actualHotspot.ref());
        assertEquals(expectedHotspot.alt(), actualHotspot.alt());
    }

    @NotNull
    private static ImmutableTransvarRecord.Builder baseRecord() {
        return ImmutableTransvarRecord.builder().transcript("irrelevant").chromosome("1").variantSpanMultipleExons(false);
    }

    @NotNull
    private static ImmutableVariantHotspotImpl.Builder baseHotspot() {
        return ImmutableVariantHotspotImpl.builder().chromosome("1");
    }
}