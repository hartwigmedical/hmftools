package com.hartwig.hmftools.common.virusbreakend;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class VirusBreakendFactoryTest {

    private static final double EPSILON = 1e-10;

    private static final String VIRUS_BREAKEND_TSV = Resources.getResource("virusbreakend/sample.virusbreakend.vcf.summary.tsv").getPath();
    private static final String VIRUS_BREAKEND_WITH_TSV = Resources.getResource(
            "virusbreakend/sample_with.virusbreakend.vcf.summary.tsv").getPath();


    @Test
    public void canReadVirusBreakendTsv() throws IOException {
        List<VirusBreakend> virusbreakendList = VirusBreakendFactory.readVirusBreakend(VIRUS_BREAKEND_TSV);

        assertEquals(0, virusbreakendList.size());

    }

    @Test
    public void canReadVirusBreakendWithTsv() throws IOException {
        List<VirusBreakend> virusbreakendList = VirusBreakendFactory.readVirusBreakend(VIRUS_BREAKEND_WITH_TSV);
        assertEquals(2, virusbreakendList.size());

        VirusBreakend virusbreakend = virusbreakendList.get(0);

        assertEquals(10, virusbreakend.taxidGenus());
        assertEquals("Alphapapillomavirus", virusbreakend.nameGenus());
        assertEquals(100, virusbreakend.readsGenusTree());
        assertEquals(110, virusbreakend.taxidSpecies());
        assertEquals("Alphapapillomavirus 9", virusbreakend.nameSpecies());
        assertEquals(11, virusbreakend.readsSpeciesTree());
        assertEquals(50, virusbreakend.taxidAssigned());
        assertEquals("Human papillomavirus type 16", virusbreakend.nameAssigned());
        assertEquals(3000, virusbreakend.readsAssignedTree());
        assertEquals(3500, virusbreakend.readsAssignedDirect());
        assertEquals("MG10.1", virusbreakend.reference());
        assertEquals(87, virusbreakend.referenceTaxid());
        assertEquals(88, virusbreakend.referenceKmerCount());
        assertEquals(1200, virusbreakend.alternateKmerCount());
        assertEquals("adjusted_MG10.1", virusbreakend.Rname());
        assertEquals(1, virusbreakend.startpos());
        assertEquals(11, virusbreakend.endpos());
        assertEquals(12, virusbreakend.numreads());
        assertEquals(100, virusbreakend.covbases());
        assertEquals(14.102, virusbreakend.coverage(), EPSILON);
        assertEquals(12.1, virusbreakend.meandepth(), EPSILON);
        assertEquals(115.2, virusbreakend.meanbaseq(), EPSILON);
        assertEquals(60, virusbreakend.meanmapq(), EPSILON);
        assertEquals(9, virusbreakend.integrations());
        assertEquals(Strings.EMPTY, virusbreakend.QCStatus());
    }
}