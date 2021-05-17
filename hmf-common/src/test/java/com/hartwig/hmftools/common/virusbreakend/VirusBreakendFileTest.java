package com.hartwig.hmftools.common.virusbreakend;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class VirusBreakendFileTest {

    private static final double EPSILON = 1e-10;

    private static final String EMPTY_VIRUS_BREAKEND_TSV =
            Resources.getResource("virusbreakend/empty.virusbreakend.vcf.summary.tsv").getPath();
    private static final String SAMPLE_VIRUS_BREAKEND_TSV =
            Resources.getResource("virusbreakend/sample.virusbreakend.vcf.summary.tsv").getPath();

    @Test
    public void canReadEmptyVirusBreakendTsv() throws IOException {
        List<VirusBreakend> virusbreakendList = VirusBreakendFile.read(EMPTY_VIRUS_BREAKEND_TSV);

        assertEquals(0, virusbreakendList.size());
    }

    @Test
    public void canReadSampleVirusBreakendTsv() throws IOException {
        List<VirusBreakend> virusbreakendList = VirusBreakendFile.read(SAMPLE_VIRUS_BREAKEND_TSV);
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
        assertEquals("adjusted_MG10.1", virusbreakend.RName());
        assertEquals(1, virusbreakend.startPos());
        assertEquals(11, virusbreakend.endPos());
        assertEquals(12, virusbreakend.numReads());
        assertEquals(100, virusbreakend.covBases());
        assertEquals(14.102, virusbreakend.coverage(), EPSILON);
        assertEquals(12.1, virusbreakend.meanDepth(), EPSILON);
        assertEquals(115.2, virusbreakend.meanBaseQ(), EPSILON);
        assertEquals(60, virusbreakend.meanMapQ(), EPSILON);
        assertEquals(9, virusbreakend.integrations());
        assertEquals(VirusBreakendQCStatus.PASS, virusbreakend.qcStatus());
    }
}