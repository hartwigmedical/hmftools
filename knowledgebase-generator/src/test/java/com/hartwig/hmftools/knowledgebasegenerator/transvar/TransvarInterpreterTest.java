package com.hartwig.hmftools.knowledgebasegenerator.transvar;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class TransvarInterpreterTest {

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

        List<VariantHotspot> hotspots = TransvarInterpreter.convertRecordToHotspots(record, Strand.REVERSE);

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