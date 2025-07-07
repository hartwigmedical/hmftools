package com.hartwig.hmftools.esvee.pon_gen;

import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.esvee.pon_gen.LocationCounter.isValid;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import org.junit.Test;

public class PonLocationsTest
{
    @Test
    public void testPonLocations()
    {
        PonLocations location = new PonLocations("id");

        location.addPosition(10);
        location.addPosition(20);
        location.addPosition(1);
        location.addPosition(5);
        location.addPosition(30);
        location.addPosition(25);
        location.addPosition(25);
        location.addPosition(2);
        location.addPosition(40);
        location.addPosition(10);
        location.addPosition(10);
        location.addPosition(5);
        location.addPosition(40);

        assertTrue(isValid(location.Locations));
        assertEquals(8, location.Locations.size());
        assertEquals(8, location.locationCount());
        assertEquals(2, location.Locations.stream().filter(x -> x.Position == 40).findFirst().orElse(null).getCount());
        assertEquals(2, location.Locations.stream().filter(x -> x.Position == 5).findFirst().orElse(null).getCount());
        assertEquals(3, location.Locations.stream().filter(x -> x.Position == 10).findFirst().orElse(null).getCount());
    }

    @Test
    public void testPonSvs()
    {
        PonStore ponStore = new PonStore();

        String chr1 = "1";
        String chr2 = "2";

        ponStore.addSv(chr1, chr1, ORIENT_FWD, ORIENT_REV, 100, 200); // DEL
        ponStore.addSv(chr1, chr1, ORIENT_FWD, ORIENT_REV, 100, 200); // DEL
        ponStore.addSv(chr1, chr1, ORIENT_FWD, ORIENT_REV, 100, 300); // another DEL

        ponStore.addSv(chr1, chr2, ORIENT_FWD, ORIENT_FWD, 100, 200); // BND

        ponStore.addSv(chr2, chr2, ORIENT_FWD, ORIENT_FWD, 100, 200); // INV

        ponStore.addSv(chr2, chr2, ORIENT_REV, ORIENT_REV, 100, 200); // another INV
        ponStore.addSv(chr2, chr2, ORIENT_REV, ORIENT_REV, 100, 200);
        ponStore.addSv(chr2, chr2, ORIENT_REV, ORIENT_REV, 100, 400); // another INV
        ponStore.addSv(chr2, chr2, ORIENT_REV, ORIENT_REV, 100, 400);
        ponStore.addSv(chr2, chr2, ORIENT_REV, ORIENT_REV, 100, 400);
        ponStore.addSv(chr2, chr2, ORIENT_REV, ORIENT_REV, 50, 200); // another INV
        ponStore.addSv(chr2, chr2, ORIENT_REV, ORIENT_REV, 50, 200);
        ponStore.addSv(chr2, chr2, ORIENT_REV, ORIENT_REV, 50, 200);
        ponStore.addSv(chr2, chr2, ORIENT_REV, ORIENT_REV, 50, 200);

        assertEquals(4, ponStore.svLocationCount());

        LocationCounter locCounter = ponStore.getLocationCounter(chr1, chr1, ORIENT_FWD, ORIENT_REV, 100, 200);
        assertTrue(locCounter != null);
        assertEquals(2, locCounter.getCount());

        locCounter = ponStore.getLocationCounter(chr1, chr1, ORIENT_FWD, ORIENT_REV, 100, 300);
        assertTrue(locCounter != null);
        assertEquals(1, locCounter.getCount());

        locCounter = ponStore.getLocationCounter(chr2, chr2, ORIENT_REV, ORIENT_REV, 100, 200);
        assertTrue(locCounter != null);
        assertEquals(2, locCounter.getCount());

        locCounter = ponStore.getLocationCounter(chr2, chr2, ORIENT_REV, ORIENT_REV, 100, 400);
        assertTrue(locCounter != null);
        assertEquals(3, locCounter.getCount());

        locCounter = ponStore.getLocationCounter(chr2, chr2, ORIENT_REV, ORIENT_REV, 50, 200);
        assertTrue(locCounter != null);
        assertEquals(4, locCounter.getCount());
    }

    @Test
    public void testPonSgls()
    {
        PonStore ponStore = new PonStore();

        String chr1 = "1";
        String chr2 = "2";

        ponStore.addSgl(chr1, ORIENT_FWD, 100);
        ponStore.addSgl(chr1, ORIENT_FWD, 100);

        ponStore.addSgl(chr1, ORIENT_REV, 100);
        ponStore.addSgl(chr1, ORIENT_REV, 100);
        ponStore.addSgl(chr1, ORIENT_REV, 100);

        ponStore.addSgl(chr1, ORIENT_FWD, 200);
        ponStore.addSgl(chr1, ORIENT_FWD, 200);
        ponStore.addSgl(chr1, ORIENT_FWD, 200);
        ponStore.addSgl(chr1, ORIENT_FWD, 200);

        ponStore.addSgl(chr2, ORIENT_FWD, 100);
        ponStore.addSgl(chr2, ORIENT_FWD, 100);
        ponStore.addSgl(chr2, ORIENT_FWD, 100);
        ponStore.addSgl(chr2, ORIENT_FWD, 100);
        ponStore.addSgl(chr2, ORIENT_FWD, 100);

        assertEquals(3, ponStore.sglLocationCount());

        LocationCounter locCounter = ponStore.getLocationCounter(chr1, ORIENT_FWD, 100);
        assertTrue(locCounter != null);
        assertEquals(2, locCounter.getCount());

        locCounter = ponStore.getLocationCounter(chr1, ORIENT_REV, 100);
        assertTrue(locCounter != null);
        assertEquals(3, locCounter.getCount());

        locCounter = ponStore.getLocationCounter(chr1, ORIENT_FWD, 200);
        assertTrue(locCounter != null);
        assertEquals(4, locCounter.getCount());

        locCounter = ponStore.getLocationCounter(chr2, ORIENT_FWD, 100);
        assertTrue(locCounter != null);
        assertEquals(5, locCounter.getCount());
    }

}
