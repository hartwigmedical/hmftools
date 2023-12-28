package com.hartwig.hmftools.esvee.processor;

import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import com.hartwig.hmftools.esvee.models.AlignedAssembly;
import com.hartwig.hmftools.esvee.models.Alignment;
import com.hartwig.hmftools.esvee.models.ExtendedAssembly;
import com.hartwig.hmftools.esvee.models.GappedAssembly;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.models.SupportedAssembly;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

public class VariantCallerTest
{
    private static final String SEQUENCE = "ATCGATCGGCTAGCTACGATCGATTACGTACGATCGATCGGCTAGCTACGATCGATTACGTACG"; // 64 bases
    private static final SAMFileHeader GERMLINE = new SAMFileHeader();
    private static final SAMFileHeader SOMATIC = new SAMFileHeader();
    private static final AtomicInteger READ_COUNTER = new AtomicInteger();

    static
    {
        GERMLINE.setAttribute("userTag", "germline");
        SOMATIC.setAttribute("userTag", "tumor");

        final SAMReadGroupRecord germlineReadGroup = new SAMReadGroupRecord("dummy");
        germlineReadGroup.setSample("germline");
        GERMLINE.addReadGroup(germlineReadGroup);

        final SAMReadGroupRecord somaticReadGroup = new SAMReadGroupRecord("dummy");
        somaticReadGroup.setSample("somatic");
        SOMATIC.addReadGroup(somaticReadGroup);
    }

    private static AlignedAssembly createAssembly(final Alignment... alignments)
    {
        return createAssembly(Arrays.asList(alignments));
    }

    private static AlignedAssembly createAssembly(final List<Alignment> alignments)
    {
        final int sequenceLength = alignments.stream()
                .mapToInt(a -> a.SequenceStartPosition + a.Length - 1)
                .max()
                .orElseThrow();
        final String bases = SEQUENCE.substring(0, sequenceLength);

        return createAssembly(bases, alignments);
    }

    private static AlignedAssembly createAssembly(final String bases, final List<Alignment> alignments)
    {
        final ExtendedAssembly dummy = new ExtendedAssembly("Dummy", bases, new SupportedAssembly("Dummy", bases));
        final var supported = new GappedAssembly("Dummy", List.of(dummy));
        supported.addEvidenceAt(createRecord(bases, GERMLINE), 0);
        supported.addEvidenceAt(createRecord(bases, SOMATIC), 0);

        return new AlignedAssembly(supported, alignments);
    }

    private static Read createRecord(final String assembly, final SAMFileHeader header)
    {
        final var record = new SAMRecord(header);
        record.setReadName("Read" + READ_COUNTER.incrementAndGet());
        record.setReadUnmappedFlag(true);
        record.setReadBases(assembly.getBytes());
        record.setBaseQualities(new byte[assembly.length()]);
        record.setAttribute("RG", "dummy");

        return new Read(record);
    }

    private static VariantCaller caller()
    {
        /*
        // TODO: Test utils
        final SVAConfig config = HMFConfig.load(Map.of(
                        "bam_file", "dummy",
                        "output_file", "dummy",
                        "ref_genome", "dummy",
                        "ref_genome_index", "dummy",
                        "junction_file", "dummy"),
                SVAConfig.class, ImmutableSVAConfig.builder());
        return new VariantCaller(config, Executors.newSingleThreadExecutor());
       */

        return null;
    }

    /* CHASHA FIXME
    @Test
    public void translocationABStandard()
    {
        // Chr 1:110 -> Chr 2:121
        // ATCGATCGGCTAGCTACGAT
        // 12345678901234567890
        //    Chr 1 || Chr 2
        //        110
        //           121

        final var assembly = createAssembly(
                new Alignment("1", 101, 1, 10, false, 60),
                new Alignment("2", 121, 11, 10, false, 60));

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).hasSize(1);
        final var variant = variants.get(0);
        assertTrue(variant.LeftChromosome).isEqualTo("1");
        assertTrue(variant.LeftPosition).isEqualTo(110);
        assertTrue(variant.RightChromosome).isEqualTo("2");
        assertTrue(variant.RightPosition).isEqualTo(121);
        assertTrue(variant.LeftDescriptor).isEqualTo("C[2:121[");
        assertTrue(variant.RightDescriptor).isEqualTo("]1:110]T");
        assertTrue(variant.Classification.toString()).isEqualTo("TRANSLOCATION");
        assertTrue(variant.quality()).isEqualTo(61);
        assertTrue(variant.germlineSupport()).isEqualTo(1);
        assertTrue(variant.somaticSupport()).isEqualTo(1);
    }

    @Test
    public void translocationABStandardInsert()
    {
        // Chr 1:110 -> Chr 2:121
        // ATCGATCGGCTAGCTACGATCGAT
        // 1234567890    1234567890
        //    Chr 1 |    | Chr 2
        //        110
        //               121

        final var assembly = createAssembly(
                new Alignment("1", 101, 1, 10, false, 60),
                Alignment.unmapped(11, 4),
                new Alignment("2", 121, 15, 10, false, 60));

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).hasSize(1);
        final var variant = variants.get(0);
        assertTrue(variant.LeftChromosome).isEqualTo("1");
        assertTrue(variant.LeftPosition).isEqualTo(110);
        assertTrue(variant.RightChromosome).isEqualTo("2");
        assertTrue(variant.RightPosition).isEqualTo(121);
        assertTrue(variant.LeftDescriptor).isEqualTo("CTAGC[2:121[");
        assertTrue(variant.RightDescriptor).isEqualTo("]1:110]TAGCT");
        assertTrue(variant.Classification.toString()).isEqualTo("TRANSLOCATION");
        assertTrue(variant.quality()).isEqualTo(61);
        assertTrue(variant.germlineSupport()).isEqualTo(1);
        assertTrue(variant.somaticSupport()).isEqualTo(1);
    }

    @Test
    public void translocationABLeftInverted()
    {
        // Chr 1:101 -> Chr 2:121
        // ATCGATCGGCTAGCTACGAT
        // 09876543211234567890
        //    Chr 1 || Chr 2
        //        101
        //           121

        final var assembly = createAssembly(
                new Alignment("1", 101, 1, 10, true, 60),
                new Alignment("2", 121, 11, 10, false, 60));

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).hasSize(1);
        final var variant = variants.get(0);
        assertTrue(variant.LeftChromosome).isEqualTo("1");
        assertTrue(variant.LeftPosition).isEqualTo(101);
        assertTrue(variant.RightChromosome).isEqualTo("2");
        assertTrue(variant.RightPosition).isEqualTo(121);
        assertTrue(variant.LeftDescriptor).isEqualTo("[2:121[G");
        assertTrue(variant.RightDescriptor).isEqualTo("[1:101[T");
        assertTrue(variant.Classification.toString()).isEqualTo("TRANSLOCATION");
        assertTrue(variant.quality()).isEqualTo(61);
        assertTrue(variant.germlineSupport()).isEqualTo(1);
        assertTrue(variant.somaticSupport()).isEqualTo(1);
    }

    @Test
    public void translocationABLeftInvertedInsert()
    {
        // Chr 1:101 -> Chr 2:121
        // ATCGATCGGCTAGCTACGATCGAT
        // 0987654321    1234567890
        //    Chr 1 |    | Chr 2
        //        101
        //               121

        final var assembly = createAssembly(
                new Alignment("1", 101, 1, 10, true, 60),
                Alignment.unmapped(11, 4),
                new Alignment("2", 121, 15, 10, false, 60));

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).hasSize(1);
        final var variant = variants.get(0);
        assertTrue(variant.LeftChromosome).isEqualTo("1");
        assertTrue(variant.LeftPosition).isEqualTo(101);
        assertTrue(variant.RightChromosome).isEqualTo("2");
        assertTrue(variant.RightPosition).isEqualTo(121);
        assertTrue(variant.LeftDescriptor).isEqualTo("[2:121[GCTAG");
        assertTrue(variant.RightDescriptor).isEqualTo("[1:101[TAGCT");
        assertTrue(variant.Classification.toString()).isEqualTo("TRANSLOCATION");
        assertTrue(variant.quality()).isEqualTo(61);
        assertTrue(variant.germlineSupport()).isEqualTo(1);
        assertTrue(variant.somaticSupport()).isEqualTo(1);
    }

    @Test
    public void translocationABRightInverted()
    {
        // Chr 1:110 -> Chr 2:130
        // ATCGATCGGCTAGCTACGAT
        // 12345678900987654321
        //    Chr 1 || Chr 2
        //        110
        //           130

        final var assembly = createAssembly(
                new Alignment("1", 101, 1, 10, false, 60),
                new Alignment("2", 121, 11, 10, true, 60));

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).hasSize(1);
        final var variant = variants.get(0);
        assertTrue(variant.LeftChromosome).isEqualTo("1");
        assertTrue(variant.LeftPosition).isEqualTo(110);
        assertTrue(variant.RightChromosome).isEqualTo("2");
        assertTrue(variant.RightPosition).isEqualTo(130);
        assertTrue(variant.LeftDescriptor).isEqualTo("C]2:130]");
        assertTrue(variant.RightDescriptor).isEqualTo("A]1:110]");
        assertTrue(variant.Classification.toString()).isEqualTo("TRANSLOCATION");
        assertTrue(variant.quality()).isEqualTo(61);
        assertTrue(variant.germlineSupport()).isEqualTo(1);
        assertTrue(variant.somaticSupport()).isEqualTo(1);
    }

    @Test
    public void translocationABRightInvertedInsert()
    {
        // Chr 1:110 -> Chr 2:130
        // ATCGATCGGCTAGCTACGATCGAT
        // 1234567890    0987654321
        //    Chr 1 |    | Chr 2
        //        110
        //               130

        final var assembly = createAssembly(
                new Alignment("1", 101, 1, 10, false, 60),
                Alignment.unmapped(11, 4),
                new Alignment("2", 121, 15, 10, true, 60));

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).hasSize(1);
        final var variant = variants.get(0);
        assertTrue(variant.LeftChromosome).isEqualTo("1");
        assertTrue(variant.LeftPosition).isEqualTo(110);
        assertTrue(variant.RightChromosome).isEqualTo("2");
        assertTrue(variant.RightPosition).isEqualTo(130);
        assertTrue(variant.LeftDescriptor).isEqualTo("CTAGC]2:130]");
        assertTrue(variant.RightDescriptor).isEqualTo("AGCTA]1:110]");
        assertTrue(variant.Classification.toString()).isEqualTo("TRANSLOCATION");
        assertTrue(variant.quality()).isEqualTo(61);
        assertTrue(variant.germlineSupport()).isEqualTo(1);
        assertTrue(variant.somaticSupport()).isEqualTo(1);
    }

    @Test
    public void translocationABBothInverted()
    {
        // Chr 1:101 -> Chr 2:130
        // ATCGATCGGCTAGCTACGAT
        // 09876543210987654321
        //    Chr 1 || Chr 2
        //        101
        //           130

        final var assembly = createAssembly(
                new Alignment("1", 101, 1, 10, true, 60),
                new Alignment("2", 121, 11, 10, true, 60));

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).hasSize(1);
        final var variant = variants.get(0);
        assertTrue(variant.LeftChromosome).isEqualTo("1");
        assertTrue(variant.LeftPosition).isEqualTo(101);
        assertTrue(variant.RightChromosome).isEqualTo("2");
        assertTrue(variant.RightPosition).isEqualTo(130);
        assertTrue(variant.LeftDescriptor).isEqualTo("]2:130]G");
        assertTrue(variant.RightDescriptor).isEqualTo("A[1:101[");
        assertTrue(variant.Classification.toString()).isEqualTo("TRANSLOCATION");
        assertTrue(variant.quality()).isEqualTo(61);
        assertTrue(variant.germlineSupport()).isEqualTo(1);
        assertTrue(variant.somaticSupport()).isEqualTo(1);
    }

    @Test
    public void translocationABBothInvertedInsert()
    {
        // Chr 1:101 -> Chr 2:130
        // ATCGATCGGCTAGCTACGATCGAT
        // 0987654321    0987654321
        //    Chr 1 |    | Chr 2
        //        101
        //               130

        final var assembly = createAssembly(
                new Alignment("1", 101, 1, 10, true, 60),
                Alignment.unmapped(11, 4),
                new Alignment("2", 121, 15, 10, true, 60));

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).hasSize(1);
        final var variant = variants.get(0);
        assertTrue(variant.LeftChromosome).isEqualTo("1");
        assertTrue(variant.LeftPosition).isEqualTo(101);
        assertTrue(variant.RightChromosome).isEqualTo("2");
        assertTrue(variant.RightPosition).isEqualTo(130);
        assertTrue(variant.LeftDescriptor).isEqualTo("]2:130]GCTAG");
        assertTrue(variant.RightDescriptor).isEqualTo("AGCTA[1:101[");
        assertTrue(variant.Classification.toString()).isEqualTo("TRANSLOCATION");
        assertTrue(variant.quality()).isEqualTo(61);
        assertTrue(variant.germlineSupport()).isEqualTo(1);
        assertTrue(variant.somaticSupport()).isEqualTo(1);
    }

    @Test
    public void translocationAALeftInverted()
    {
        // Chr 1:101 -> Chr 1:121
        // ATCGATCGGCTAGCTACGAT
        // 09876543211234567890
        //    Chr 1 || Chr 1
        //        101
        //           121

        final var assembly = createAssembly(
                new Alignment("1", 101, 1, 10, true, 60),
                new Alignment("1", 121, 11, 10, false, 60));

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).hasSize(1);
        final var variant = variants.get(0);
        assertTrue(variant.LeftChromosome).isEqualTo("1");
        assertTrue(variant.LeftPosition).isEqualTo(101);
        assertTrue(variant.RightChromosome).isEqualTo("1");
        assertTrue(variant.RightPosition).isEqualTo(121);
        assertTrue(variant.LeftDescriptor).isEqualTo("[1:121[G");
        assertTrue(variant.RightDescriptor).isEqualTo("[1:101[T");
        assertTrue(variant.Classification.toString()).isEqualTo("INVERSION");
        assertTrue(variant.quality()).isEqualTo(61);
        assertTrue(variant.germlineSupport()).isEqualTo(1);
        assertTrue(variant.somaticSupport()).isEqualTo(1);
    }

    @Test
    public void translocationAALeftInvertedInsert()
    {
        // Chr 1:101 -> Chr 1:121
        // ATCGATCGGCTAGCTACGATCGAT
        // 0987654321    1234567890
        //    Chr 1 |    | Chr 1
        //        101
        //               121

        final var assembly = createAssembly(
                new Alignment("1", 101, 1, 10, true, 60),
                Alignment.unmapped(11, 4),
                new Alignment("1", 121, 15, 10, false, 60));

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).hasSize(1);
        final var variant = variants.get(0);
        assertTrue(variant.LeftChromosome).isEqualTo("1");
        assertTrue(variant.LeftPosition).isEqualTo(101);
        assertTrue(variant.RightChromosome).isEqualTo("1");
        assertTrue(variant.RightPosition).isEqualTo(121);
        assertTrue(variant.LeftDescriptor).isEqualTo("[1:121[GCTAG");
        assertTrue(variant.RightDescriptor).isEqualTo("[1:101[TAGCT");
        assertTrue(variant.Classification.toString()).isEqualTo("INVERSION");
        assertTrue(variant.quality()).isEqualTo(61);
        assertTrue(variant.germlineSupport()).isEqualTo(1);
        assertTrue(variant.somaticSupport()).isEqualTo(1);
    }

    @Test
    public void translocationAARightInverted()
    {
        // Chr 1:110 -> Chr 1:130
        // ATCGATCGGCTAGCTACGAT
        // 12345678900987654321
        //    Chr 1 || Chr 1
        //        110
        //           130

        final var assembly = createAssembly(
                new Alignment("1", 101, 1, 10, false, 60),
                new Alignment("1", 121, 11, 10, true, 60));

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).hasSize(1);
        final var variant = variants.get(0);
        assertTrue(variant.LeftChromosome).isEqualTo("1");
        assertTrue(variant.LeftPosition).isEqualTo(110);
        assertTrue(variant.RightChromosome).isEqualTo("1");
        assertTrue(variant.RightPosition).isEqualTo(130);
        assertTrue(variant.LeftDescriptor).isEqualTo("C]1:130]");
        assertTrue(variant.RightDescriptor).isEqualTo("A]1:110]");
        assertTrue(variant.Classification.toString()).isEqualTo("INVERSION");
        assertTrue(variant.quality()).isEqualTo(61);
        assertTrue(variant.germlineSupport()).isEqualTo(1);
        assertTrue(variant.somaticSupport()).isEqualTo(1);
    }

    @Test
    public void translocationAARightInvertedInsert()
    {
        // Chr 1:110 -> Chr 1:130
        // ATCGATCGGCTAGCTACGATCGAT
        // 1234567890    0987654321
        //    Chr 1 |    | Chr 1
        //        110
        //               130

        final var assembly = createAssembly(
                new Alignment("1", 101, 1, 10, false, 60),
                Alignment.unmapped(11, 4),
                new Alignment("1", 121, 15, 10, true, 60));

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).hasSize(1);
        final var variant = variants.get(0);
        assertTrue(variant.LeftChromosome).isEqualTo("1");
        assertTrue(variant.LeftPosition).isEqualTo(110);
        assertTrue(variant.RightChromosome).isEqualTo("1");
        assertTrue(variant.RightPosition).isEqualTo(130);
        assertTrue(variant.LeftDescriptor).isEqualTo("CTAGC]1:130]");
        assertTrue(variant.RightDescriptor).isEqualTo("AGCTA]1:110]");
        assertTrue(variant.Classification.toString()).isEqualTo("INVERSION");
        assertTrue(variant.quality()).isEqualTo(61);
        assertTrue(variant.germlineSupport()).isEqualTo(1);
        assertTrue(variant.somaticSupport()).isEqualTo(1);
    }

    @Test
    public void deleteShort()
    {
        // Chr 1:110 -> Chr 1:121
        // ATCGATCGGCTAGCTACGAT
        // 12345678901234567890
        //    Chr 1 || Chr 1

        final var assembly = createAssembly(
                new Alignment("1", 101, 1, 10, false, 60),
                new Alignment("1", 121, 11, 10, false, 60));

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).isEmpty();
    }

    @Test
    public void deleteStandard()
    {
        // Chr 1:110 -> Chr 1:221
        // ATCGATCGGCTAGCTACGAT
        // 12345678901234567890
        //    Chr 1 || Chr 1

        final var assembly = createAssembly(
                new Alignment("1", 101, 1, 10, false, 60),
                new Alignment("1", 221, 11, 10, false, 60));

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).hasSize(1);
        final var variant = variants.get(0);
        assertTrue(variant.LeftChromosome).isEqualTo("1");
        assertTrue(variant.LeftPosition).isEqualTo(110);
        assertTrue(variant.RightChromosome).isEqualTo("1");
        assertTrue(variant.RightPosition).isEqualTo(221);
        assertTrue(variant.LeftDescriptor).isEqualTo("C[1:221[");
        assertTrue(variant.RightDescriptor).isEqualTo("]1:110]T");
        assertTrue(variant.Classification.toString()).isEqualTo("DELETION 110");
        assertTrue(variant.quality()).isEqualTo(61);
        assertTrue(variant.germlineSupport()).isEqualTo(1);
        assertTrue(variant.somaticSupport()).isEqualTo(1);
    }

    @Test
    public void deleteStandardInsert()
    {
        // Chr 1:110 -> Chr 1:221
        // ATCGATCGGCTAGCTACGATCGAT
        // 1234567890    1234567890
        //    Chr 1 |    | Chr 1

        final var assembly = createAssembly(
                new Alignment("1", 101, 1, 10, false, 60),
                Alignment.unmapped(11, 4),
                new Alignment("1", 221, 15, 10, false, 60));

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).hasSize(1);
        final var variant = variants.get(0);
        assertTrue(variant.LeftChromosome).isEqualTo("1");
        assertTrue(variant.LeftPosition).isEqualTo(110);
        assertTrue(variant.RightChromosome).isEqualTo("1");
        assertTrue(variant.RightPosition).isEqualTo(221);
        assertTrue(variant.LeftDescriptor).isEqualTo("CTAGC[1:221[");
        assertTrue(variant.RightDescriptor).isEqualTo("]1:110]TAGCT");
        assertTrue(variant.Classification.toString()).isEqualTo("DELETION 110");
        assertTrue(variant.quality()).isEqualTo(61);
        assertTrue(variant.germlineSupport()).isEqualTo(1);
        assertTrue(variant.somaticSupport()).isEqualTo(1);
    }

    @Test
    public void insertExact()
    {
        // Chr 1:110 -> Chr 1:111
        // ATCGATCGGCTAGCTACGATCGATTACGTACGATCGATCGGCTAGCTACGATCGATTACG
        //           1234567980123456798012345679801234567980
        // 1234567890                                        1234567890
        //    Chr 1 |                                        | Chr 1

        final var assembly = createAssembly(
                new Alignment("1", 101, 1, 10, false, 60),
                Alignment.insert(11, 40),
                new Alignment("1", 111, 51, 10, false, 60));

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).hasSize(1);
        final var variant = variants.get(0);
        assertTrue(variant.LeftChromosome).isEqualTo("1");
        assertTrue(variant.LeftPosition).isEqualTo(110);
        assertTrue(variant.RightChromosome).isEqualTo("1");
        assertTrue(variant.RightPosition).isEqualTo(111);
        assertTrue(variant.LeftDescriptor).isEqualTo("CTAGCTACGATCGATTACGTACGATCGATCGGCTAGCTACG[1:111[");
        assertTrue(variant.RightDescriptor).isEqualTo("]1:110]TAGCTACGATCGATTACGTACGATCGATCGGCTAGCTACGA");
        assertTrue(variant.Classification.toString()).isEqualTo("INSERT 40");
        assertTrue(variant.quality()).isEqualTo(61);
        assertTrue(variant.germlineSupport()).isEqualTo(1);
        assertTrue(variant.somaticSupport()).isEqualTo(1);
    }

    @Test
    public void insertSmallGap()
    {
        // Chr 1:110 -> Chr 1:113
        // ATCGATCGGCTAGCTACGATCGATTACGTACGATCGATCGGCTAGCTACGATCGATTACG
        //           1234567980123456798012345679801234567980
        // 1234567890                                        3456789012
        //    Chr 1 |                                        | Chr 1

        final var assembly = createAssembly(
                new Alignment("1", 101, 1, 10, false, 60),
                Alignment.insert(11, 40),
                new Alignment("1", 113, 51, 10, false, 60));

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).hasSize(1);
        final var variant = variants.get(0);
        assertTrue(variant.LeftChromosome).isEqualTo("1");
        assertTrue(variant.LeftPosition).isEqualTo(110);
        assertTrue(variant.RightChromosome).isEqualTo("1");
        assertTrue(variant.RightPosition).isEqualTo(113);
        assertTrue(variant.LeftDescriptor).isEqualTo("CTAGCTACGATCGATTACGTACGATCGATCGGCTAGCTACG[1:113[");
        assertTrue(variant.RightDescriptor).isEqualTo("]1:110]TAGCTACGATCGATTACGTACGATCGATCGGCTAGCTACGA");
        assertTrue(variant.Classification.toString()).isEqualTo("INSERT 40");
        assertTrue(variant.quality()).isEqualTo(61);
        assertTrue(variant.germlineSupport()).isEqualTo(1);
        assertTrue(variant.somaticSupport()).isEqualTo(1);
    }

    @Test
    public void insertSmallOverlap()
    {
        // Chr 1:110 -> Chr 1:109
        // ATCGATCGGCTAGCTACGATCGATTACGTACGATCGATCGGCTAGCTACGATCGATTACG
        //           1234567980123456798012345679801234567980
        // 1234567890                                        9012345678
        //    Chr 1 |                                        | Chr 1

        final var assembly = createAssembly(
                new Alignment("1", 101, 1, 10, false, 60),
                Alignment.insert(11, 40),
                new Alignment("1", 109, 51, 10, false, 60));

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).hasSize(1);
        final var variant = variants.get(0);
        assertTrue(variant.LeftChromosome).isEqualTo("1");
        assertTrue(variant.LeftPosition).isEqualTo(109);
        assertTrue(variant.RightChromosome).isEqualTo("1");
        assertTrue(variant.RightPosition).isEqualTo(110);
        assertTrue(variant.LeftDescriptor).isEqualTo("]1:110]TAGCTACGATCGATTACGTACGATCGATCGGCTAGCTACGA");
        assertTrue(variant.RightDescriptor).isEqualTo("CTAGCTACGATCGATTACGTACGATCGATCGGCTAGCTACG[1:109[");
        assertTrue(variant.Classification.toString()).isEqualTo("INSERT 40");
        assertTrue(variant.quality()).isEqualTo(61);
        assertTrue(variant.germlineSupport()).isEqualTo(1);
        assertTrue(variant.somaticSupport()).isEqualTo(1);
    }

    @Test
    public void duplicationNoInsert()
    {
        // Chr 1:110 -> Chr 1:50
        // ATCGATCGGCTAGCTACGAT
        // 12345678900123456789
        //    Chr 1 || Chr 1

        final var assembly = createAssembly(
                new Alignment("1", 101, 1, 10, false, 60),
                new Alignment("1", 50, 11, 10, false, 60));

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).hasSize(1);
        final var variant = variants.get(0);
        assertTrue(variant.LeftChromosome).isEqualTo("1");
        assertTrue(variant.LeftPosition).isEqualTo(50);
        assertTrue(variant.RightChromosome).isEqualTo("1");
        assertTrue(variant.RightPosition).isEqualTo(110);
        assertTrue(variant.LeftDescriptor).isEqualTo("]1:110]T");
        assertTrue(variant.RightDescriptor).isEqualTo("C[1:50[");
        assertTrue(variant.Classification.toString()).isEqualTo("DUPLICATION 61");
        assertTrue(variant.quality()).isEqualTo(61);
        assertTrue(variant.germlineSupport()).isEqualTo(1);
        assertTrue(variant.somaticSupport()).isEqualTo(1);
    }

    @Test
    public void callsOverhangCorrectly()
    {
        // Chr 1:110 -> Chr 1:111
        // ATCGATCGGCTAGCTACGATCGATTACGTACGATCGATCGGCTAGCTACGATCGATTACG
        //           1234567980123456798012345679801234567980
        // 1234567890                                        1234567890
        //    Chr 1 |                                        | Chr 1

        final var assembly = new AlignedAssembly(new GappedAssembly("dummy",
                List.of(new ExtendedAssembly("dummy", SEQUENCE, null))),
                List.of(
                        new Alignment("1", 101, 1, 10, false, 60),
                        Alignment.insert(11, 40),
                        new Alignment("1", 111, 51, 10, false, 60)
                ));
        assembly.addEvidenceAt(createRecord("CGATCGGCTAGCTACGATCGATTACGTACGATCGATCGGCTAGCTACGATCGAT", SOMATIC), 2);

        final var variants = caller().callVariants(List.of(assembly));
        assertTrue(variants).hasSize(1);
        final var variant = variants.get(0);
        assertTrue(variant.LeftChromosome).isEqualTo("1");
        assertTrue(variant.LeftPosition).isEqualTo(110);
        assertTrue(variant.RightChromosome).isEqualTo("1");
        assertTrue(variant.RightPosition).isEqualTo(111);
        assertTrue(variant.LeftDescriptor).isEqualTo("CTAGCTACGATCGATTACGTACGATCGATCGGCTAGCTACG[1:111[");
        assertTrue(variant.RightDescriptor).isEqualTo("]1:110]TAGCTACGATCGATTACGTACGATCGATCGGCTAGCTACGA");
        assertTrue(variant.Classification.toString()).isEqualTo("INSERT 40");
        assertTrue(variant.germlineSupport()).isEqualTo(0);
        assertTrue(variant.somaticSupport()).isEqualTo(1);
        assertTrue(variant.variantAssemblies()).hasSize(1);
        final VariantCall.VariantAssembly variantAssembly = variant.variantAssemblies().iterator().next();
        assertTrue(variantAssembly.LeftOverhang).isEqualTo(8);
        assertTrue(variantAssembly.RightOverhang).isEqualTo(6);
    }
    */
}
