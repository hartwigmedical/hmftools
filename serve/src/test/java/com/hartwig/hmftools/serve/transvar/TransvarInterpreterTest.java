package com.hartwig.hmftools.serve.transvar;

import static com.hartwig.hmftools.serve.transvar.TransvarTestFactory.testInterpreter;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarComplexInsertDelete;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarDeletion;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarDuplication;
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
                .annotation(ImmutableTransvarSnvMnv.builder()
                        .gdnaRef("TA")
                        .gdnaAlt("GC")
                        .referenceCodon("TAC")
                        .addCandidateCodons("GCA", "GCC", "GCG", "GCT")
                        .build())
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
                .annotation(ImmutableTransvarSnvMnv.builder()
                        .gdnaRef("CA")
                        .gdnaAlt("GC")
                        .referenceCodon("TGG")
                        .addCandidateCodons("GCA", "GCC", "GCG", "GCT")
                        .build())
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
        TransvarRecord record =
                baseRecord().gdnaPosition(5).annotation(ImmutableTransvarDeletion.builder().deletedBases("GAT").build()).build();

        List<VariantHotspot> hotspots = testInterpreter().convertRecordToHotspots(record, Strand.FORWARD);

        assertEquals(1, hotspots.size());

        assertHotspot(baseHotspot().position(4).ref("CGAT").alt("C").build(), hotspots.get(0));
    }

    @Test
    public void canConvertInsertionToHotspots() {
        TransvarRecord forwardRecord =
                baseRecord().gdnaPosition(5).annotation(ImmutableTransvarInsertion.builder().insertedBases("GAA").build()).build();

        List<VariantHotspot> forwardHotspots = testInterpreter().convertRecordToHotspots(forwardRecord, Strand.FORWARD);

        assertEquals(2, forwardHotspots.size());

        assertHotspot(baseHotspot().position(4).ref("C").alt("CGAG").build(), forwardHotspots.get(0));
        assertHotspot(baseHotspot().position(4).ref("C").alt("CGAA").build(), forwardHotspots.get(1));

        TransvarRecord reverseRecord =
                baseRecord().gdnaPosition(5).annotation(ImmutableTransvarInsertion.builder().insertedBases("GAA").build()).build();

        List<VariantHotspot> reverseHotspots = testInterpreter().convertRecordToHotspots(reverseRecord, Strand.REVERSE);

        assertEquals(2, reverseHotspots.size());

        assertHotspot(baseHotspot().position(4).ref("C").alt("CAAA").build(), reverseHotspots.get(0));
        assertHotspot(baseHotspot().position(4).ref("C").alt("CGAA").build(), reverseHotspots.get(1));
    }

    @Test
    public void canConvertComplexDeletionInsertionToHotspot() {
        TransvarRecord complexRecord = baseRecord().gdnaPosition(2)
                .annotation(ImmutableTransvarComplexInsertDelete.builder()
                        .deletedBaseCount(4)
                        .addCandidateAlternativeSequences("GGG")
                        .addCandidateAlternativeSequences("CCC")
                        .build())
                .build();

        List<VariantHotspot> hotspots = testInterpreter().convertRecordToHotspots(complexRecord, Strand.FORWARD);

        assertEquals(2, hotspots.size());

        assertHotspot(baseHotspot().position(1).ref("GATCG").alt("GGGG").build(), hotspots.get(0));
        assertHotspot(baseHotspot().position(1).ref("GATCG").alt("GCCC").build(), hotspots.get(1));
    }

    @Test
    public void canConvertDuplicationToHotspot() {
        TransvarRecord record =
                baseRecord().gdnaPosition(5).annotation(ImmutableTransvarDuplication.builder().duplicatedBaseCount(3).build()).build();

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
}