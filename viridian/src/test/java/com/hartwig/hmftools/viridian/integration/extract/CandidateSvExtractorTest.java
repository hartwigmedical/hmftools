package com.hartwig.hmftools.viridian.integration.extract;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.INTEGRATION_INSERT_LENGTH_MIN;

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

public class CandidateSvExtractorTest
{
    private static final String TUMOR_ID = "TUMOR";

    // Inserted sequences at and just below the single length threshold, plus a longer one.
    private static final String INSERT_AT_MIN_LENGTH = "ACGT".repeat(20).substring(0, INTEGRATION_INSERT_LENGTH_MIN);
    private static final String INSERT_BELOW_MIN_LENGTH = INSERT_AT_MIN_LENGTH.substring(1);
    private static final String LONG_INSERT = "ACGT".repeat(15);

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

    @Test
    public void testExtractSingleBreakendWithAnnotations()
    {
        // A single breakend carrying every annotation the phase records but never acts on.

        String records = record(
                "chr1", 1000, "sgl_1", "A", "A" + INSERT_AT_MIN_LENGTH + ".", "minQual",
                "SVTYPE=SGL;LINE;INSALN=chr7:100|+|50M|60;INSRMRC=SINE;INSRMRT=Alu;INSRMP=0.8",
                "0:30:10:0.0", "12:40:20:0.25");

        BreakendSupport support = new BreakendSupport(12, 60, 0.25, 0, 40);
        CandidateSv expected = new CandidateSv(
                StructuralVariantType.SGL, "minQual",
                new HostBreakend("sgl_1", new BasePosition("chr1", 1000), FORWARD, support), null,
                INSERT_AT_MIN_LENGTH, true, new InsertRepeat("SINE", "Alu", 0.8), "chr7:100|+|50M|60");

        assertEquals(List.of(expected), extract(records));
    }

    @Test
    public void testExtractPairedSvAsOneCandidate()
    {
        // A paired variant yields one candidate carrying both host breakends, not one per breakend. It is identified by the
        // id ESVEE shares between the two records, rather than by either record's own id.

        String records = pairedDeletion("del_1", LONG_INSERT);

        BreakendSupport startSupport = new BreakendSupport(8, 55, 0.15, 0, 45);
        BreakendSupport endSupport = new BreakendSupport(8, 57, 0.15, 0, 47);
        CandidateSv expected = new CandidateSv(
                StructuralVariantType.DEL, "PASS",
                new HostBreakend("del_1_o", new BasePosition("chr1", 5000), FORWARD, startSupport),
                new HostBreakend("del_1_h", new BasePosition("chr1", 9000), REVERSE, endSupport),
                LONG_INSERT, false, null, "");

        assertEquals(List.of(expected), extract(records));
    }

    @Test
    public void testExtractAppliesInsertLengthThreshold()
    {
        // Insert length is the only candidate gate, and the floor is the same for every variant type: a paired variant at the
        // floor qualifies just as a single breakend does.

        String records = singleBreakend("chr1", 1000, "sgl_short", INSERT_BELOW_MIN_LENGTH)
                + singleBreakend("chr1", 2000, "sgl_long", INSERT_AT_MIN_LENGTH)
                + pairedDeletion("del_short", INSERT_BELOW_MIN_LENGTH)
                + pairedDeletion("del_long", INSERT_AT_MIN_LENGTH);

        assertEquals(List.of("sgl_long", "del_long_o"), breakendIds(extract(records)));
    }

    @Test
    public void testExtractSkipsNonHumanContigs()
    {
        // Viral sequence never becomes a breakend coordinate, so a non-human contig is host-irrelevant and skipped.

        String records = singleBreakend("chrEBV", 1000, "ebv_1", INSERT_AT_MIN_LENGTH);

        assertEquals(List.of(), extract(records));
    }

    @Test
    public void testExtractSkipsBreakendWithoutMate()
    {
        String records = record(
                "chr1", 5000, "orphan_o", "A", "A" + LONG_INSERT + "[chr1:9000[", "PASS",
                "SVTYPE=DEL;MATEID=never_arrives", "0:1:1:0.0", "1:1:1:0.1");

        assertEquals(List.of(), extract(records));
    }

    @Test
    public void testExtractResolvesTumorGenotypeByName()
    {
        // Fragment counts must follow the named tumor sample, not a genotype ordinal, so the columns are ordered
        // tumor-then-reference here to catch an ordinal assumption.

        String header = HEADER.replace("NORMAL\tTUMOR", "TUMOR\tNORMAL");
        String records = record(
                "chr1", 1000, "sgl_1", "A", "A" + INSERT_AT_MIN_LENGTH + ".", "PASS", "SVTYPE=SGL", "12:40:20:0.25", "0:30:10:0.0");

        List<CandidateSv> candidates = extract(header, records);
        assertEquals(new BreakendSupport(12, 60, 0.25, 0, 40), candidates.get(0).startBreakend().support());
    }

    @Test
    public void testExtractFailsWhenTumorSampleAbsent()
    {
        String header = HEADER.replace("NORMAL\tTUMOR", "NORMAL\tOTHER");
        String records = singleBreakend("chr1", 1000, "sgl_1", INSERT_AT_MIN_LENGTH);

        assertThrows(UserInputError.class, () -> extract(header, records));
    }

    private static String singleBreakend(String chromosome, int position, String id, String insertSequence)
    {
        return record(
                chromosome, position, id, "A", "A" + insertSequence + ".", "PASS", "SVTYPE=SGL", "0:1:1:0.0", "1:1:1:0.1");
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

    private static List<String> breakendIds(List<CandidateSv> candidates)
    {
        return candidates.stream().map(candidate -> candidate.startBreakend().id()).toList();
    }

    private List<CandidateSv> extract(String records)
    {
        return extract(HEADER, records);
    }

    private List<CandidateSv> extract(String header, String records)
    {
        try
        {
            File vcfFile = mTempDir.newFile("esvee.unfiltered.vcf");
            Files.writeString(vcfFile.toPath(), header + records);
            return new CandidateSvExtractor(TUMOR_ID).extract(vcfFile.getAbsolutePath());
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
    }
}
