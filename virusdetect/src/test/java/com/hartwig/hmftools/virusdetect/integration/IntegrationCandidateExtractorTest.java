package com.hartwig.hmftools.virusdetect.integration;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.virusdetect.UserInputError;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class IntegrationCandidateExtractorTest
{
    private static final String TUMOR_ID = "TUMOR";

    // Inserted sequences at and either side of the two length thresholds.
    private static final String INSERT_20 = "ACGTACGTACGTACGTACGT";
    private static final String INSERT_19 = "ACGTACGTACGTACGTACG";
    private static final String INSERT_50 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
    private static final String INSERT_49 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA";

    // Tumor genotype is second here, matching the usual reference-then-tumor ordering.
    private static final String HEADER = """
            ##fileformat=VCFv4.2
            ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="">
            ##INFO=<ID=MATEID,Number=1,Type=String,Description="">
            ##INFO=<ID=LINE,Number=0,Type=Flag,Description="">
            ##INFO=<ID=INSALN,Number=1,Type=String,Description="">
            ##INFO=<ID=INSRMRC,Number=1,Type=String,Description="">
            ##INFO=<ID=INSRMRT,Number=1,Type=String,Description="">
            ##INFO=<ID=INSRMRO,Number=1,Type=String,Description="">
            ##INFO=<ID=INSRMP,Number=1,Type=Float,Description="">
            ##FORMAT=<ID=VF,Number=1,Type=Integer,Description="">
            ##FORMAT=<ID=REF,Number=1,Type=Integer,Description="">
            ##FORMAT=<ID=REFPAIR,Number=1,Type=Integer,Description="">
            ##FORMAT=<ID=AF,Number=1,Type=Float,Description="">
            ##FILTER=<ID=minQual,Description="">
            ##contig=<ID=chr1,length=248956422>
            ##contig=<ID=chr2,length=242193529>
            ##contig=<ID=chrEBV,length=171823>
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR
            """;

    private static final String FORMAT = "VF:REF:REFPAIR:AF";

    @Rule
    public TemporaryFolder mTempDir = new TemporaryFolder();

    // A single breakend carrying every annotation the phase records but never acts on.
    @Test
    public void extractsSingleBreakendWithAnnotations()
    {
        String records = record(
                "chr1", 1000, "sgl_1", "A", "A" + INSERT_20 + ".", "minQual",
                "SVTYPE=SGL;LINE;INSALN=chr7:100|+|50M|60;INSRMRC=SINE;INSRMRT=Alu;INSRMRO=-;INSRMP=0.8",
                "0:30:10:0.0", "12:40:20:0.25");

        BreakendSupport support = new BreakendSupport(12, 60, 0.25, 0, 40);
        IntegrationCandidate expected = new IntegrationCandidate(
                "sgl_1", StructuralVariantType.SGL, "minQual",
                new HostBreakend("chr1", 1000, FORWARD, support), null,
                INSERT_20, true, "SINE", "Alu", (byte) -1, 0.8, "chr7:100|+|50M|60");

        assertEquals(List.of(expected), extract(records));
    }

    // A paired variant yields one candidate carrying both host breakends, not one per breakend.
    // It is identified by its start breakend's record id.
    @Test
    public void extractsPairedVariantAsOneCandidate()
    {
        String records = pairedDeletion("del_1", INSERT_50);

        BreakendSupport startSupport = new BreakendSupport(8, 55, 0.15, 0, 45);
        BreakendSupport endSupport = new BreakendSupport(8, 57, 0.15, 0, 47);
        IntegrationCandidate expected = new IntegrationCandidate(
                "del_1_o", StructuralVariantType.DEL, "PASS",
                new HostBreakend("chr1", 5000, FORWARD, startSupport),
                new HostBreakend("chr1", 9000, REVERSE, endSupport),
                INSERT_50, false, null, null, null, null, "");

        assertEquals(List.of(expected), extract(records));
    }

    // Insert length is the only gate: a single breakend needs 20 bases, a paired variant 50.
    @Test
    public void appliesInsertLengthThresholds()
    {
        String records = record(
                "chr1", 1000, "sgl_short", "A", "A" + INSERT_19 + ".", "PASS", "SVTYPE=SGL", "0:1:1:0.0", "1:1:1:0.1")
                + record(
                "chr1", 2000, "sgl_long", "A", "A" + INSERT_20 + ".", "PASS", "SVTYPE=SGL", "0:1:1:0.0", "1:1:1:0.1")
                + pairedDeletion("del_short", INSERT_49)
                + pairedDeletion("del_long", INSERT_50);

        assertEquals(List.of("sgl_long", "del_long_o"), svIds(extract(records)));
    }

    // Viral sequence never becomes a breakend coordinate, so a non-human contig is host-irrelevant and skipped.
    @Test
    public void skipsNonHumanContigs()
    {
        String records = record(
                "chrEBV", 1000, "ebv_1", "A", "A" + INSERT_20 + ".", "PASS", "SVTYPE=SGL", "0:1:1:0.0", "1:1:1:0.1");

        assertEquals(List.of(), extract(records));
    }

    // A paired breakend whose mate never arrives is incomplete, so it cannot become a candidate.
    @Test
    public void skipsBreakendWithoutMate()
    {
        String records = record(
                "chr1", 5000, "orphan_o", "A", "A" + INSERT_50 + "[chr1:9000[", "PASS",
                "SVTYPE=DEL;MATEID=never_arrives", "0:1:1:0.0", "1:1:1:0.1");

        assertEquals(List.of(), extract(records));
    }

    // Fragment counts must follow the named tumor sample, not a genotype ordinal, so the columns are ordered
    // tumor-then-reference here to catch an ordinal assumption.
    @Test
    public void resolvesTumorGenotypeByName()
    {
        String header = HEADER.replace("NORMAL\tTUMOR", "TUMOR\tNORMAL");
        String records = record(
                "chr1", 1000, "sgl_1", "A", "A" + INSERT_20 + ".", "PASS", "SVTYPE=SGL", "12:40:20:0.25", "0:30:10:0.0");

        List<IntegrationCandidate> candidates = extract(header, records);
        assertEquals(new BreakendSupport(12, 60, 0.25, 0, 40), candidates.get(0).startBreakend().support());
    }

    @Test
    public void failsWhenTumorSampleAbsent()
    {
        String header = HEADER.replace("NORMAL\tTUMOR", "NORMAL\tOTHER");
        String records = record(
                "chr1", 1000, "sgl_1", "A", "A" + INSERT_20 + ".", "PASS", "SVTYPE=SGL", "0:1:1:0.0", "1:1:1:0.1");

        assertThrows(UserInputError.class, () -> extract(header, records));
    }

    private static String pairedDeletion(String id, String insertSequence)
    {
        return record(
                "chr1", 5000, id + "_o", "A", "A" + insertSequence + "[chr1:9000[", "PASS",
                "SVTYPE=DEL;MATEID=" + id + "_h", "0:45:0:0.0", "8:55:0:0.15")
                + record(
                "chr1", 9000, id + "_h", "A", "]chr1:5000]" + insertSequence + "A", "PASS",
                "SVTYPE=DEL;MATEID=" + id + "_o", "0:47:0:0.0", "8:57:0:0.15");
    }

    private static String record(
            String chromosome, int position, String id, String ref, String alt, String filter, String info,
            String normalGenotype, String tumorGenotype)
    {
        return String.join(
                "\t", chromosome, String.valueOf(position), id, ref, alt, "100", filter, info, FORMAT,
                normalGenotype, tumorGenotype) + "\n";
    }

    private static List<String> svIds(List<IntegrationCandidate> candidates)
    {
        return candidates.stream().map(IntegrationCandidate::svId).toList();
    }

    private List<IntegrationCandidate> extract(String records)
    {
        return extract(HEADER, records);
    }

    private List<IntegrationCandidate> extract(String header, String records)
    {
        try
        {
            File vcfFile = mTempDir.newFile("esvee.unfiltered.vcf");
            Files.writeString(vcfFile.toPath(), header + records);
            return new IntegrationCandidateExtractor(TUMOR_ID).extract(vcfFile.getAbsolutePath());
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
    }
}
