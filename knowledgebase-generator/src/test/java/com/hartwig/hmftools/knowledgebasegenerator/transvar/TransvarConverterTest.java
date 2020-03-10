package com.hartwig.hmftools.knowledgebasegenerator.transvar;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import org.junit.Test;

public class TransvarConverterTest {

    @Test
    public void canConvertSnvToRecord() {
        String line =
                "MTOR:p.L2230V\tENST00000361445 (protein_coding)\tMTOR\t-\tchr1:g.11182158A>C/c.6688T>G/p.L2230V\tinside_[cds_in_exon_48]"
                        + "\tCSQN=Missense;reference_codon=TTA;candidate_codons=GTA,GTC,GTG,GTT;candidate_mnv_variants="
                        + "chr1:g.11182156_11182158delTAAinsGAC,chr1:g.11182156_11182158delTAAinsCAC,chr1:g.11182156_11182158delTAAinsAAC;"
                        + "aliases=ENSP00000354558;source=Ensembl";

        TransvarRecord record = TransvarConverter.toTransvarRecord(line);

        assertNotNull(record);
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

    @Test
    public void canConvertMnvToRecord() {
        String line = "TET2:p.Y1294A\tENST00000540549 (protein_coding)\tTET2\t+\t"
                + "chr4:g.106180852_106180853delTAinsGC/c.3880_3881delTAinsGC/p.Y1294A\tinside_[cds_in_exon_7]\t"
                + "CSQN=Missense;reference_codon=TAC;candidate_codons=GCA,GCC,GCG,GCT;candidate_mnv_variants="
                + "chr4:g.106180852_106180854delTACinsGCA,chr4:g.106180852_106180854delTACinsGCG,chr4:"
                + "g.106180852_106180854delTACinsGCT;aliases=ENSP00000442788;source=Ensembl";

        TransvarRecord record = TransvarConverter.toTransvarRecord(line);

        assertNotNull(record);
        assertEquals("ENST00000540549", record.transcript());
        assertEquals("4", record.chromosome());
        assertEquals(106180852, record.gdnaPosition());
        assertEquals("TA", record.gdnaRef());
        assertEquals("GC", record.gdnaAlt());
        assertEquals("TAC", record.referenceCodon());
        assertEquals("GCA", record.candidateCodons().get(0));
        assertEquals("GCC", record.candidateCodons().get(1));
        assertEquals("GCG", record.candidateCodons().get(2));
        assertEquals("GCT", record.candidateCodons().get(3));
    }

    @Test
    public void canConvertDeletionToRecord() {
        String deletionLine =
                "NOTCH1:p.V1578del\tENST00000277541 (protein_coding)\tNOTCH1\t-\tchr9:g.139399420_139399422delCCA/c.4732_4734delGTG/"
                        + "p.V1578delV\tinside_[cds_in_exon_26]\tCSQN=InFrameDeletion;left_align_gDNA=g.139399409_139399411delCAC;"
                        + "unaligned_gDNA=g.139399409_139399411delCAC;left_align_cDNA=c.4721_4723delTGG;unalign_cDNA=c.4732_4734delGTG;"
                        + "left_align_protein=p.V1575delV;unalign_protein=p.V1578delV;imprecise;aliases=ENSP00000277541;source=Ensembl";

        TransvarRecord deletion = TransvarConverter.toTransvarRecord(deletionLine);

        assertNotNull(deletion);
        assertEquals("ENST00000277541", deletion.transcript());
        assertEquals("9", deletion.chromosome());
        assertEquals(139399420, deletion.gdnaPosition());
        assertEquals("CCA", deletion.gdnaRef());
        assertEquals("", deletion.gdnaAlt());
    }

    @Test
    public void canConvertInsertionToRecord() {
        String insertionLine =
                "ERBB2:p.G776_V777insYVMA\tENST00000584450 (protein_coding)\tERBB2\t+\tchr17:g.37880999_37881000insTATGTAATGGCA/"
                        + "c.2328_2329insTATGTAATGGCA/p.G776_V777insYVMA\tinside_[cds_in_exon_20]\tCSQN=InFrameInsertion;"
                        + "left_align_protein=p.G776_V777insYVMA;unalign_protein=p.G776_V777insYVMA;left_align_gDNA="
                        + "g.37880999_37881000insTATGTAATGGCA;unalign_gDNA=g.37880999_37881000insTATGTAATGGCA;left_align_cDNA="
                        + "c.2328_2329insTATGTAATGGCA;unalign_cDNA=c.2328_2329insTATGTAATGGCA;32_CandidatesOmitted;"
                        + "aliases=ENSP00000463714;source=Ensembl";

        TransvarRecord insertion = TransvarConverter.toTransvarRecord(insertionLine);

        assertNotNull(insertion);
        assertEquals("ENST00000584450", insertion.transcript());
        assertEquals("17", insertion.chromosome());
        assertEquals(37880999, insertion.gdnaPosition());
        assertEquals("", insertion.gdnaRef());
        assertEquals("TATGTAATGGCA", insertion.gdnaAlt());
    }

    @Test
    public void canConvertDuplicationToRecord() {
        String dupLine = "ERBB2:p.Y772_A775dup\tENST00000584450 (protein_coding)\tERBB2\t+\tchr17:g.37880985_37880996/"
                + "c.2314_2325/p.Y772_A775\tinside_[cds_in_exon_20]\tprotein_sequence=YVMA;cDNA_sequence=TAC..GCT;"
                + "gDNA_sequence=TAC..GCT;aliases=ENSP00000463714;source=Ensembl";

        TransvarRecord duplication = TransvarConverter.toTransvarRecord(dupLine);

        assertNotNull(duplication);
        assertEquals("ENST00000584450", duplication.transcript());
        assertEquals("17", duplication.chromosome());
        assertEquals(37880985, duplication.gdnaPosition());
        assertEquals("", duplication.gdnaRef());
        assertEquals("", duplication.gdnaAlt());
        assertEquals(12, (int) duplication.dupLength());
    }

    @Test
    public void unresolvedProteinAnnotationLeadsToNull() {
        String line = "RAD50:p.L1273F\t.\t.\t.\t././.\t.\tno_valid_transcript_found";

        assertNull(TransvarConverter.toTransvarRecord(line));
    }
}