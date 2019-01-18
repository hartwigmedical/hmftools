package com.hartwig.hmftools.svanalysis.visualisation;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.junit.Test;

public class SvRegionFileTest {

    private static final String REGIONS_PATH = Resources.getResource("visualisation/ClusterLinks.tsv").getPath();
    private static final String LINKS_PATH = Resources.getResource("visualisation/ClusterSV.tsv").getPath();


    @Test
    public void testCOLO829() throws IOException {
        final String prefix = "/Users/jon/hmf/analysis/sv/SvWriter";

        final List<GenomeRegion> regions = SvRegionFile.readLinks(REGIONS_PATH);
        final List<SvLink> links = SvLinkFile.readLinks(LINKS_PATH);


        SvCircosWriter writer = new SvCircosWriter(prefix);
        writer.writeLinks(regions, links);
    }

}
