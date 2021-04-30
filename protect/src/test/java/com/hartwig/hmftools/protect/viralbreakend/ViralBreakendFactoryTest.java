package com.hartwig.hmftools.protect.viralbreakend;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.protect.viralbreakend.ViralBreakendFactory;
import com.hartwig.hmftools.protect.viralbreakend.Viralbreakend;

import org.junit.Test;

public class ViralBreakendFactoryTest {

    private static final String VIRAL_BREAKEND_TSV = Resources.getResource("test_run/viralbreakend/sample.virusbreakend.vcf.summary.tsv").getPath();
    private static final String VIRAL_BREAKEND_WITH_TSV = Resources.getResource(
            "test_run/viralbreakend/sample_with.virusbreakend.vcf.summary.tsv").getPath();


    @Test
    public void canReadViralBreakendTsv() throws IOException {
        List<Viralbreakend> viralbreakendList = ViralBreakendFactory.readViralBreakend(VIRAL_BREAKEND_TSV);

        assertEquals(0, viralbreakendList.size());

    }

    @Test
    public void canReadViralBreakendWithTsv() throws IOException {
        List<Viralbreakend> viralbreakendList = ViralBreakendFactory.readViralBreakend(VIRAL_BREAKEND_WITH_TSV);
        assertEquals(2, viralbreakendList.size());

        Viralbreakend viralbreakend = viralbreakendList.get(0);

        assertEquals("10", viralbreakend.taxidGenus());
        assertEquals("Alphapapillomavirus", viralbreakend.nameGenus());
        assertEquals("100", viralbreakend.readsGenusTree());
        assertEquals("1200", viralbreakend.taxidSpecies());
        assertEquals("Alphapapillomavirus 9", viralbreakend.nameSpecies());
        assertEquals("5", viralbreakend.readsSpeciesTree());
        assertEquals("450", viralbreakend.taxidAssigned());
        assertEquals("Human papillomavirus type 16", viralbreakend.nameAssigned());
        assertEquals("12", viralbreakend.readsAssignedTree());
        assertEquals("47", viralbreakend.readsAssignedDirect());
        assertEquals("AB1234.1", viralbreakend.Reference());
        assertEquals("3567", viralbreakend.referenceTaxid());
        assertEquals("768", viralbreakend.referenceKmerCount());
        assertEquals("54321", viralbreakend.alternateKmerCountRname());
        assertEquals("adjusted_AB1234.1", viralbreakend.startpos());
        assertEquals("12", viralbreakend.endpos());
        assertEquals("89", viralbreakend.numreads());
        assertEquals("123", viralbreakend.covbases());
        assertEquals("789", viralbreakend.coverage());
        assertEquals("1.1", viralbreakend.meandepth());
        assertEquals("2.1", viralbreakend.meanbaseq());
        assertEquals("3.2", viralbreakend.meanmapq());
        assertEquals("16", viralbreakend.integrations());
        assertEquals("3", viralbreakend.QCStatus());
    }
}