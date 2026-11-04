package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.common.sv.SvVcfTags.INSALN;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.esvee.caller.CallerTestUtils.createSv;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.assembly.alignment.AlternativeAlignment;

import org.junit.Test;

public class HotspotTest
{
    @Test
    public void testHotspotClassification()
    {
        HotspotCache hotspotCache = new HotspotCache(null);

        ChrBaseRegion regionStart = new ChrBaseRegion(CHR_1, 1000, 10000);
        ChrBaseRegion regionEnd = new ChrBaseRegion(CHR_2, 20000, 30000);
        hotspotCache.addHotspot(new KnownHotspot(regionStart, Orientation.FORWARD, regionEnd, Orientation.REVERSE, "KNOWN_PAIR"));

        String vcfId = "ID_01";
        String insSeq = "";

        Variant var = createSv(
                vcfId, CHR_1, CHR_2, 1050, 25000, Orientation.FORWARD.asByte(), Orientation.REVERSE.asByte(), insSeq);

        assertTrue(hotspotCache.isHotspotVariant(var));

        // invalid orientation
        var = createSv(vcfId, CHR_1, CHR_2, 1050, 25000, Orientation.REVERSE.asByte(), Orientation.REVERSE.asByte(), insSeq);

        assertFalse(hotspotCache.isHotspotVariant(var));

        // SGLs
        Map<String,Object> commonAttributes = Maps.newHashMap();

        var = createSv(
                vcfId, CHR_1, null, 1050, 0, ORIENT_FWD, ORIENT_REV, "",
                commonAttributes, Collections.emptyMap(), Collections.emptyMap());

        assertFalse(hotspotCache.isHotspotVariant(var));

        AlternativeAlignment altAlignment = new AlternativeAlignment(CHR_2, 25000, Orientation.REVERSE, "100M", 0);

        commonAttributes.put(INSALN, altAlignment.vcfString());

        var = createSv(
                vcfId, CHR_1, null, 1050, 0, ORIENT_FWD, ORIENT_REV, "",
                commonAttributes, Collections.emptyMap(), Collections.emptyMap());

        assertTrue(hotspotCache.isHotspotVariant(var));
    }
}
