package com.hartwig.hmftools.serve.transvar;

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
        assertNull(record.indelLength());
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
        assertNull(record.indelLength());
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
        assertNull(deletion.referenceCodon());
        assertNull(deletion.candidateCodons());
        assertEquals(3, (int) deletion.indelLength());
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
        assertNull(insertion.referenceCodon());
        assertNull(insertion.candidateCodons());
        assertEquals(12, (int) insertion.indelLength());
    }

    @Test
    public void canConvertInsertionAndDeletionToRecord() {
        // TODO
        String delInsLine = "EGFR:p.L747_A750delinsP\tENST00000275493 (protein_coding)\tEGFR\t+\tchr7:g.55242469_55242480delinsCCT/"
                + "c.2239_2250delinsCCT/p.L747_A750delinsP\tinside_[cds_in_exon_19]\tCSQN=MultiAAMissense;" +
                "candidate_alternative_sequence=CCT/CCG/CCA/CCC;aliases=ENSP00000275493;source=Ensembl";

        TransvarRecord delins = TransvarConverter.toTransvarRecord(delInsLine);

        assertNotNull(delins);
        assertEquals("ENST00000275493", delins.transcript());
        assertEquals("7", delins.chromosome());
        assertEquals(55242469, delins.gdnaPosition());
    }

    @Test
    public void canConvertDuplicationToRecord() {
        String dupLineNoBases = "ERBB2:p.Y772_A775dup\tENST00000584450 (protein_coding)\tERBB2\t+\tchr17:g.37880985_37880996/"
                + "c.2314_2325/p.Y772_A775\tinside_[cds_in_exon_20]\tprotein_sequence=YVMA;cDNA_sequence=TAC..GCT;"
                + "gDNA_sequence=TAC..GCT;aliases=ENSP00000463714;source=Ensembl";

        TransvarRecord duplication1 = TransvarConverter.toTransvarRecord(dupLineNoBases);

        assertNotNull(duplication1);
        assertEquals("ENST00000584450", duplication1.transcript());
        assertEquals("17", duplication1.chromosome());
        assertEquals(37880985, duplication1.gdnaPosition());
        assertEquals("", duplication1.gdnaRef());
        assertEquals("", duplication1.gdnaAlt());
        assertNull(duplication1.referenceCodon());
        assertNull(duplication1.candidateCodons());
        assertEquals(12, (int) duplication1.indelLength());

        String dupLineWithBases = "BRAF:p.T599_V600insV\tENST00000288602 (protein_coding)\tBRAF\t-\tchr7:g.140453136_140453138dupACT"
                + "/c.1797_1799dupAGT/p.V600dupV\tinside_[cds_in_exon_15]\tCSQN=InFrameInsertion;left_align_protein=p.T599_V600insV;"
                + "unalign_protein=p.T599_V600insV;left_align_gDNA=g.140453135_140453136insACT;unalign_gDNA=g.140453137_140453138insTAC;"
                + "left_align_cDNA=c.1796_1797insAGT;unalign_cDNA=c.1797_1798insGTA;4_CandidatesOmitted;aliases=ENSP00000288602;"
                + "source=Ensembl";

        TransvarRecord duplication2 = TransvarConverter.toTransvarRecord(dupLineWithBases);

        assertNotNull(duplication2);
        assertEquals("ENST00000288602", duplication2.transcript());
        assertEquals("7", duplication2.chromosome());
        assertEquals(140453136, duplication2.gdnaPosition());
        assertEquals("", duplication2.gdnaRef());
        assertEquals("", duplication2.gdnaAlt());
        assertNull(duplication2.referenceCodon());
        assertNull(duplication2.candidateCodons());
        assertEquals(3, (int) duplication2.indelLength());
    }

    @Test
    public void longRangeLeadsToNull() {
        String line = "BRCA2:p.V1839_E1901del\tENST00000544455 (protein_coding)\tBRCA2\t+\tchr13:g.32914008_32914196del189/"
                + "c.5516_5704del189/p.V1839_E1901del63\tinside_[cds_in_exon_11]\tCSQN=InFrameDeletion;left_align_gDNA="
                + "g.32914004_32914192del189;unaligned_gDNA=g.32914007_32914195del189;left_align_cDNA=c.5512_5700del189;"
                + "unalign_cDNA=c.5515_5703del189;left_align_protein=p.E1838_S1900del63;unalign_protein=p.V1839_E1901del63;"
                + "imprecise;aliases=ENSP00000439902;source=Ensembl";

        assertNull(TransvarConverter.toTransvarRecord(line));
    }

    @Test
    public void unresolvedProteinAnnotationLeadsToNull() {
        String line = "RAD50:p.L1273F\t.\t.\t.\t././.\t.\tno_valid_transcript_found";

        assertNull(TransvarConverter.toTransvarRecord(line));
    }

    @Test
    public void invalidMutationStringLeadsToNull() {
        String line = "ERBB2:p.Exon 20 insertions\t.\t.\t.\t././.\t.\tError_invalid_mutation_string_p.Exon 20 insertions";

        assertNull(TransvarConverter.toTransvarRecord(line));
    }

    @Test
    public void noneTypeLeadsToNull() {
        String line = "BRAF:p.T599insTT\t.\t.\t.\t././.\t.\tError_int() argument must be a string or a number, not 'NoneType'";

        assertNull(TransvarConverter.toTransvarRecord(line));
    }
}