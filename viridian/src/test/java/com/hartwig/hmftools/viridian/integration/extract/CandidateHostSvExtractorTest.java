package com.hartwig.hmftools.viridian.integration.extract;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.viridian.common.UserInputError;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class CandidateHostSvExtractorTest
{
    private static final String TUMOR_ID = "TUMOR";

    // Inserted sequences at and just below the single length threshold, plus a longer one.
    private static final String INSERT_30 = "ACGTACGTACGTACGTACGTACGTACGTAC";
    private static final String INSERT_29 = "ACGTACGTACGTACGTACGTACGTACGTA";
    private static final String INSERT_50 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";

    // Tumor genotype is second here, matching the usual reference-then-tumor ordering.
    private static final String HEADER = """
            ##fileformat=VCFv4.2
            ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="">
            ##INFO=<ID=MATEID,Number=1,Type=String,Description="">
            ##INFO=<ID=SVID,Number=1,Type=String,Description="">
            ##INFO=<ID=LINE,Number=0,Type=Flag,Description="">
            ##INFO=<ID=INSALN,Number=1,Type=String,Description="">
            ##INFO=<ID=INSRMRC,Number=1,Type=String,Description="">
            ##INFO=<ID=INSRMRT,Number=1,Type=String,Description="">
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
    public void testExtractSingleBreakendWithAnnotations()
    {
        String records = record(
                "chr1", 1000, "sgl_1", "A", "A" + INSERT_30 + ".", "minQual",
                "SVTYPE=SGL;LINE;INSALN=chr7:100|+|50M|60;INSRMRC=SINE;INSRMRT=Alu;INSRMP=0.8",
                "0:30:10:0.0", "12:40:20:0.25");

        BreakendSupport support = new BreakendSupport(12, 60, 0.25, 0, 40);
        CandidateHostSv expected = new CandidateHostSv(
                StructuralVariantType.SGL, "minQual",
                new HostBreakend("sgl_1", new BasePosition("chr1", 1000), FORWARD, support), null,
                INSERT_30, true, new InsertRepeat("SINE", "Alu", 0.8), "chr7:100|+|50M|60");

        assertEquals(List.of(expected), extract(records));
    }

    // A paired variant yields one candidate carrying both host breakends, not one per breakend. It is identified by the
    // id ESVEE shares between the two records, rather than by either record's own id.
    @Test
    public void testExtractPairedSvAsOneCandidate()
    {
        String records = pairedDeletion("del_1", INSERT_50);

        BreakendSupport startSupport = new BreakendSupport(8, 55, 0.15, 0, 45);
        BreakendSupport endSupport = new BreakendSupport(8, 57, 0.15, 0, 47);
        CandidateHostSv expected = new CandidateHostSv(
                StructuralVariantType.DEL, "PASS",
                new HostBreakend("del_1_o", new BasePosition("chr1", 5000), FORWARD, startSupport),
                new HostBreakend("del_1_h", new BasePosition("chr1", 9000), REVERSE, endSupport),
                INSERT_50, false, null, "");

        assertEquals(List.of(expected), extract(records));
    }

    // Insert length is the only candidate gate, and the floor is the same for every variant type: a paired variant at the
    // floor qualifies just as a single breakend does.
    @Test
    public void testExtractAppliesInsertLengthThreshold()
    {
        String records = record(
                "chr1", 1000, "sgl_short", "A", "A" + INSERT_29 + ".", "PASS", "SVTYPE=SGL", "0:1:1:0.0", "1:1:1:0.1")
                + record(
                "chr1", 2000, "sgl_long", "A", "A" + INSERT_30 + ".", "PASS", "SVTYPE=SGL", "0:1:1:0.0", "1:1:1:0.1")
                + pairedDeletion("del_short", INSERT_29)
                + pairedDeletion("del_long", INSERT_30);

        assertEquals(List.of("sgl_long", "del_long_o"), breakendIds(extract(records)));
    }

    // Viral sequence never becomes a breakend coordinate, so a non-human contig is host-irrelevant and skipped.
    @Test
    public void testExtractSkipsNonHumanContigs()
    {
        String records = record(
                "chrEBV", 1000, "ebv_1", "A", "A" + INSERT_30 + ".", "PASS", "SVTYPE=SGL", "0:1:1:0.0", "1:1:1:0.1");

        assertEquals(List.of(), extract(records));
    }

    // A paired breakend whose mate never arrives is incomplete, so it cannot become a candidate.
    @Test
    public void testExtractSkipsBreakendWithoutMate()
    {
        String records = record(
                "chr1", 5000, "orphan_o", "A", "A" + INSERT_50 + "[chr1:9000[", "PASS",
                "SVTYPE=DEL;MATEID=never_arrives", "0:1:1:0.0", "1:1:1:0.1");

        assertEquals(List.of(), extract(records));
    }

    // Fragment counts must follow the named tumor sample, not a genotype ordinal, so the columns are ordered
    // tumor-then-reference here to catch an ordinal assumption.
    @Test
    public void testExtractResolvesTumorGenotypeByName()
    {
        String header = HEADER.replace("NORMAL\tTUMOR", "TUMOR\tNORMAL");
        String records = record(
                "chr1", 1000, "sgl_1", "A", "A" + INSERT_30 + ".", "PASS", "SVTYPE=SGL", "12:40:20:0.25", "0:30:10:0.0");

        List<CandidateHostSv> candidates = extract(header, records);
        assertEquals(new BreakendSupport(12, 60, 0.25, 0, 40), candidates.get(0).startBreakend().support());
    }

    @Test
    public void testExtractFailsWhenTumorSampleAbsent()
    {
        String header = HEADER.replace("NORMAL\tTUMOR", "NORMAL\tOTHER");
        String records = record(
                "chr1", 1000, "sgl_1", "A", "A" + INSERT_30 + ".", "PASS", "SVTYPE=SGL", "0:1:1:0.0", "1:1:1:0.1");

        assertThrows(UserInputError.class, () -> extract(header, records));
    }

    private static String pairedDeletion(String id, String insertSequence)
    {
        return record(
                "chr1", 5000, id + "_o", "A", "A" + insertSequence + "[chr1:9000[", "PASS",
                "SVTYPE=DEL;MATEID=" + id + "_h;SVID=" + id, "0:45:0:0.0", "8:55:0:0.15")
                + record(
                "chr1", 9000, id + "_h", "A", "]chr1:5000]" + insertSequence + "A", "PASS",
                "SVTYPE=DEL;MATEID=" + id + "_o;SVID=" + id, "0:47:0:0.0", "8:57:0:0.15");
    }

    private static String record(
            String chromosome, int position, String id, String ref, String alt, String filter, String info,
            String normalGenotype, String tumorGenotype)
    {
        return String.join(
                "\t", chromosome, String.valueOf(position), id, ref, alt, "100", filter, info, FORMAT,
                normalGenotype, tumorGenotype) + "\n";
    }

    private static List<String> breakendIds(List<CandidateHostSv> candidates)
    {
        return candidates.stream().map(candidate -> candidate.startBreakend().id()).toList();
    }

    private List<CandidateHostSv> extract(String records)
    {
        return extract(HEADER, records);
    }

    private List<CandidateHostSv> extract(String header, String records)
    {
        try
        {
            File vcfFile = mTempDir.newFile("esvee.unfiltered.vcf");
            Files.writeString(vcfFile.toPath(), header + records);
            return new CandidateHostSvExtractor(TUMOR_ID).extract(vcfFile.getAbsolutePath());
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
    }
}
