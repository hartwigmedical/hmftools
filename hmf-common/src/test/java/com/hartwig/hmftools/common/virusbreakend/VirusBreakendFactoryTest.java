package com.hartwig.hmftools.common.virusbreakend;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class VirusBreakendFactoryTest {

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
        assertEquals(1200, virusbreakend.taxidSpecies());
        assertEquals("Alphapapillomavirus 9", virusbreakend.nameSpecies());
        assertEquals(5, virusbreakend.readsSpeciesTree());
        assertEquals(450, virusbreakend.taxidAssigned());
        assertEquals("Human papillomavirus type 16", virusbreakend.nameAssigned());
        assertEquals(12, virusbreakend.readsAssignedTree());
        assertEquals(47, virusbreakend.readsAssignedDirect());
        assertEquals("AB1234.1", virusbreakend.reference());
        assertEquals(3567, virusbreakend.referenceTaxid());
        assertEquals(768, virusbreakend.referenceKmerCount());
        assertEquals("adjusted_AB1234.1", virusbreakend.alternateKmerCountRname());
        assertEquals(54321, virusbreakend.startpos());
        assertEquals(12, virusbreakend.endpos());
        assertEquals(89, virusbreakend.numreads());
        assertEquals(123, virusbreakend.covbases());
        assertEquals(789, virusbreakend.coverage());
        assertEquals(1, virusbreakend.meandepth());
        assertEquals(2, virusbreakend.meanbaseq());
        assertEquals(3, virusbreakend.meanmapq());
        assertEquals(16, virusbreakend.integrations());
        assertEquals(3, virusbreakend.QCStatus());
    }
}