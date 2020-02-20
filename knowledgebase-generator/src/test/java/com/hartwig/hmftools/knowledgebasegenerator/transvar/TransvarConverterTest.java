package com.hartwig.hmftools.knowledgebasegenerator.transvar;

import static org.junit.Assert.assertEquals;

import org.junit.Ignore;
import org.junit.Test;

public class TransvarConverterTest {

    @Test
    @Ignore
    public void canConvertTransvar() {
        String line =
                "MTOR:p.L2230V\tENST00000361445 (protein_coding)\tMTOR\t-\tchr1:g.11182158A>C/c.6688T>G/p.L2230V\tinside_[cds_in_exon_48]"
                        + "\tCSQN=Missense;reference_codon=TTA;candidate_codons=GTA,GTC,GTG,GTT;candidate_mnv_variants=" +
                        "chr1:g.11182156_11182158delTAAinsGAC,chr1:g.11182156_11182158delTAAinsCAC,chr1:g.11182156_11182158delTAAinsAAC;" +
                        "aliases=ENSP00000354558;source=Ensembl";

        TransvarRecord record = TransvarConverter.toTransvarRecord(line);

        assertEquals("ENST00000361445", record.transcript());
        assertEquals("1", record.chromosome());
        assertEquals(11182158, record.gdnaPosition());
        assertEquals("A", record.gdnaRef());
        assertEquals("C", record.gdnaAlt());
        assertEquals("TTA", record.referenceCodon());
        assertEquals("GTA", record.candidateCodons().get(0));
        assertEquals("GTC", record.candidateCodons().get(1));
        assertEquals("GTG", record.candidateCodons().get(2));
        assertEquals("GTT", record.candidateCodons().get(3));
    }
}