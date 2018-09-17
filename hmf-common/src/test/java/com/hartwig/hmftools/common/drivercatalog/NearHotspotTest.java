package com.hartwig.hmftools.common.drivercatalog;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

public class NearHotspotTest {

    private NearHotspot victim;

    @Before
    public void setup() {
        victim = new NearHotspot(5, Lists.newArrayList(create("1", 1000)));
    }

    @Test
    public void testSNP() {
        assertEquals(false, victim.isNearHotspot(snp(994)));

        assertEquals(true, victim.isNearHotspot(snp(995)));
        assertEquals(true, victim.isNearHotspot(snp(996)));
        assertEquals(true, victim.isNearHotspot(snp(997)));
        assertEquals(true, victim.isNearHotspot(snp(998)));
        assertEquals(true, victim.isNearHotspot(snp(999)));

        assertEquals(true, victim.isNearHotspot(snp(1000)));

        assertEquals(true, victim.isNearHotspot(snp(1001)));
        assertEquals(true, victim.isNearHotspot(snp(1002)));
        assertEquals(true, victim.isNearHotspot(snp(1003)));
        assertEquals(true, victim.isNearHotspot(snp(1004)));
        assertEquals(true, victim.isNearHotspot(snp(1005)));

        assertEquals(false, victim.isNearHotspot(snp(1006)));
    }

    @Test
    public void testDel() {
        assertEquals(false, victim.isNearHotspot(del(993)));

        assertEquals(true, victim.isNearHotspot(del(994)));
        assertEquals(true, victim.isNearHotspot(del(995)));
        assertEquals(true, victim.isNearHotspot(del(996)));
        assertEquals(true, victim.isNearHotspot(del(997)));
        assertEquals(true, victim.isNearHotspot(del(998)));
        assertEquals(true, victim.isNearHotspot(del(999)));

        assertEquals(true, victim.isNearHotspot(snp(1000)));

        assertEquals(true, victim.isNearHotspot(del(1001)));
        assertEquals(true, victim.isNearHotspot(del(1002)));
        assertEquals(true, victim.isNearHotspot(del(1003)));
        assertEquals(true, victim.isNearHotspot(del(1004)));
        assertEquals(true, victim.isNearHotspot(del(1005)));

        assertEquals(false, victim.isNearHotspot(del(1006)));
    }

    private GenomePosition create(String chromosome, long position) {
        return new GenomePosition() {
            @NotNull
            @Override
            public String chromosome() {
                return chromosome;
            }

            @Override
            public long position() {
                return position;
            }
        };
    }

    @NotNull
    private SomaticVariant snp(long position) {
        return SomaticVariantTestBuilderFactory.create().chromosome("1").position(position).ref("A").build();
    }

    @NotNull
    private SomaticVariant del(long position) {
        return SomaticVariantTestBuilderFactory.create().chromosome("1").position(position).ref("AA").build();
    }

}
