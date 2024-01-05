package com.hartwig.hmftools.esvee.old;

import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.esvee.common.Junction;

import org.jetbrains.annotations.Nullable;
import org.junit.Ignore;

@Ignore
public class SVAssemblyIntegrationTest
{
    public static final String SOMATIC_BAM = "/Users/james/code/data/COLO829/sv-prep2/COLO829v003T.sv_prep_sorted.bam";
    public static final String GERMLINE_BAM = "/Users/james/code/data/COLO829/sv-prep2/COLO829v003R.sv_prep_sorted.bam";
    public static final String JUNCTION_FILE = "/Users/james/code/data/COLO829/sv-prep2/COLO829v003T.sv_prep.junctions.tsv";
    public static final String REFERENCE_GENOME = "/Users/james/code/data/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta";
    public static final String REFERENCE_GENOME_INDEX = "/Users/james/code/data/refgenome/Homo_sapiens.GRCh37.GATK.illumina.img";

    private int containsNearJunction(final Map<String, Set<Integer>> positions, final Junction junction, final int tolerance)
    {
        return containsNearJunction(positions, junction.chromosome(), junction.position(), tolerance);
    }

    private int containsNearJunction(final Map<String, Set<Integer>> positions, final String chromosome, final int position,
            final int tolerance)
    {
        @Nullable
        final Set<Integer> positionsInChromosome = positions.get(chromosome);
        if(positionsInChromosome == null)
            return -1;

        for(int i = 0; i <= tolerance; i++)
            if(positionsInChromosome.contains(position + i))
                return position + i;
            else if(positionsInChromosome.contains(position - i))
                return position - i;

        return -1;
    }

    /*
    public Set<VariantCall> process(final List<Pair<String, Integer>> locations, final boolean debug)
    {
        final Map<String, Object> options = new HashMap<>(Map.of(
                "bam_file", SOMATIC_BAM,
                "germline_bam", GERMLINE_BAM,
                "output_file", "output.vcf",
                "ref_genome", REFERENCE_GENOME,
                "ref_genome_index", REFERENCE_GENOME_INDEX,
                "junction_file", JUNCTION_FILE
        ));
        if(debug)
        {
            options.put("threads", "1");
            options.put("create_html_summaries", "true");
            options.put("debug", "true");
            options.put("timeouts_enabled", "false");
        }

        final SVAConfig config = HMFConfig.load(options, SVAConfig.class, ImmutableSVAConfig.builder());

        // assertTrue(config != null);

        final Context context = Context.create(config);
        final Processor processor = new Processor(context);

        final Map<String, Set<Integer>> toProcessPositions = new HashMap<>();
        for(var pair : locations)
            toProcessPositions.computeIfAbsent(pair.getKey(), ignored -> new HashSet<>()).add(pair.getValue());

        final List<Junction> junctions = JunctionReader.readJunctionFile(config.junctionFile())
                .filter(junction -> containsNearJunction(toProcessPositions, junction, 5) != -1)
                .collect(Collectors.toList());

        return new HashSet<>(processor.run(junctions));
    }
    */

    /* CHASHA FIXME
    @Test
    public void junction3_6_6_3()
    {
        final List<Pair<String, Integer>> locations = List.of(
                Pair.of("6", 26194039),
                Pair.of("6", 26194118),
                Pair.of("6", 26431917),
                Pair.of("6", 26194406),
                Pair.of("3", 26431917)
        );

        final Set<VariantCall> results = process(locations, false);

        final Set<ExpectedVariant> expected = Set.of(
                ExpectedVariant.of("6", 26194117, "G]6:26194406]", "6", 26194406, "A]6:26194117]"),
                ExpectedVariant.of("6", 26194039, ".ATCGGGTCACTCCCACCCGAATATTGCGCTTTTCAGACCGGCTTAAAAAACGGCGCACCACGAGACTATATCCC"),
                ExpectedVariant.of("3", 26431917, "[6:26194041[T", "6", 26194041, "[3:26431917[A")
        );

        assertVariants(results, expected);
    }

    @Test
    public void junction3_12_10()
    {
        final List<Pair<String, Integer>> locations = List.of(
                Pair.of("3", 25400287),
                Pair.of("3", 25400457),
                Pair.of("3", 25400602),
                Pair.of("3", 25400860),
                Pair.of("3", 25401059),
                Pair.of("12", 72666892),
                Pair.of("12", 72667074),
                Pair.of("10", 60477224),
                Pair.of("10", 60477422)
        );

        final Set<VariantCall> results = process(locations, false);

        final Set<ExpectedVariant> expected = Set.of(
                ExpectedVariant.of("3", 25401059, "G[10:60477224[", "10", 60477224, "]3:25401059]T"),
                ExpectedVariant.of("3", 25400602, "GTGAATCCATCA[12:72666892[", "12", 72666892, "]3:25400602]TGAATCCATCAG"),
                ExpectedVariant.of("10", 60477422, "GC]12:72667074]", "12", 72667074, "GG]10:60477422]")
        );

        assertVariants(results, expected);
    }

    @Ignore // Non-deterministic :(
    @Test
    public void chromosome7()
    {
        final List<Pair<String, Integer>> locations = List.of(
                Pair.of("7", 125746123),
                Pair.of("7", 126098487),
                Pair.of("7", 126166900),
                Pair.of("7", 126167440),
                Pair.of("7", 126167444),
                Pair.of("7", 143936532),
                Pair.of("7", 143936535),
                Pair.of("7", 143937207),
                Pair.of("7", 143938341),
                Pair.of("7", 143939546),
                Pair.of("7", 143939814),
                Pair.of("7", 143959226),
                Pair.of("7", 144088794),
                Pair.of("7", 144090465)
        );

        final Set<VariantCall> results = process(locations, false);

        final Set<ExpectedVariant> expected = Set.of(
                ExpectedVariant.of("7", 126167444, "GGCTTCTCTTTAGGTGTTCCAAAGGCCAGAGACTATGTCTGTAACTATTATACTAAGTTATACCAATCTCCTCCAGCCTGGACATATGGAGTCTTTTCTACCCATGATTAAAAGCTTGAGAAAAATGAGCCTGAATTTTTGACTGGCCTGAAATTTCATCCTAAGAAAGAGGATAGTAATGCAGCTAATATATTTTAACATAGGTGAAAAAACCAAAAGGGTAATACACAAATCCATGAGAAAAATGGACAAAAATATGTAGATACAGACAACAGGAATCCTCCTATGTCCAGTAGAAGAAAGAAGTAGAAATAAACATGAAAATACCTTTGACCAGCAAATAT."),
                ExpectedVariant.of("7", 126098488, "[7:126167444[T", "7", 126167444, "[7:126098488[T"),
                ExpectedVariant.of("7", 125746123, "AT[7:126166901[", "7", 126166901, "]7:125746123]TG"),
                ExpectedVariant.of("7", 144090465, "AATATATTTTACAAAATAAAACTTTACATCTTGCTCATGTAAAGTCCTAGGCAGGGATGACAGACAGCCTTCCACTGGTGATGCAGGGACCCAGGCACCTTCCATCGCATGGTCCCGCCATCCTAGAGCCCTGTAGGCCTTCACTTCCAGCCAGTGATGGGAAGGCCAAACATGGAGGAGGCATGCCCTCCTCCTCCATGCATGCCCTCCTTCCACTGCAGCCCTGGAAGCTGCACA.")
        );

        assertVariants(results, expected);
    }

    @Test
    public void simpleDelete()
    {
        final List<Pair<String, Integer>> locations = List.of(
                Pair.of("20", 1648840),
                Pair.of("20", 1648847),
                Pair.of("20", 1648864),
                Pair.of("20", 1648881)
        );

        final Set<VariantCall> results = process(locations, false);

        final Set<ExpectedVariant> expected = Set.of(
                ExpectedVariant.of("20", 1648847, "T[20:1648881[", "20", 1648881, "]20:1648847]C")
        );
        assertVariants(results, expected);
    }

    @Test
    public void simpleDelete2()
    {
        final List<Pair<String, Integer>> locations = List.of(
                Pair.of("18", 56214208),
                Pair.of("18", 56214241),
                Pair.of("18", 56214277),
                Pair.of("18", 56214279)
        );

        final Set<VariantCall> results = process(locations, false);

        final Set<ExpectedVariant> expected = Set.of(
                ExpectedVariant.of("18", 56214244, "G[18:56214280[", "18", 56214280, "]18:56214244]T")
        );
        assertVariants(results, expected);
    }

    @Test
    public void deleteShowsNoFalseDupe()
    {
        final List<Pair<String, Integer>> locations = List.of(
                Pair.of("20", 14962959),
                Pair.of("20", 14988825),
                Pair.of("20", 	15000624),
                Pair.of("20", 15013947),
                Pair.of("20", 15013842),
                Pair.of("20", 15108817)
        );

        final Set<VariantCall> results = process(locations, false);

        final Set<ExpectedVariant> expected = Set.of(
                ExpectedVariant.of("20", 15000624, "T[20:15013844[", "20", 15013844, "]20:15000624]G"),
                ExpectedVariant.of("20", 14962958, "C[20:15013951[", "20", 15013951, "]20:14962958]C"),
                ExpectedVariant.of("20", 14988825, "T[20:15108819[", "20", 15108819, "]20:14988825]A")
        );
        assertVariants(results, expected);
    }

    // This test fails if we don't correctly recognise support like:
    // ATCGGCTATCAAAAAAAAAAAAAAAACAATCGATCGA
    // ATCGGCTATCAAAAAAAAAAAAAAAAAAATCGATCGA
    // <p>
    // If we fail to see the latter as supporting the former, we will have FPs in this area.
    //
    @Ignore
    @Test
    public void ch14_inversionNearRepeat()
    {
        final List<Pair<String, Integer>> locations = List.of(
                Pair.of("14", 87902497),
                Pair.of("14", 87902499),
                Pair.of("14", 87902538)
        );

        final Set<VariantCall> results = process(locations, false);

        final Set<ExpectedVariant> expected = Set.of(
                ExpectedVariant.of("14", 87902500, "[14:87902538[T", "14", 87902538, "[14:87902500[T"),
                ExpectedVariant.of("14", 87902499, "A]14:87902538]", "14", 87902538, "T]14:87902499]")
        );
        assertVariants(results, expected);
    }

    @Test
    public void chrX_deletionBadAssembly()
    {
        final List<Pair<String, Integer>> locations = List.of(
                Pair.of("X", 34062705),
                Pair.of("X", 34059778)
        );

        final Set<VariantCall> results = process(locations, true);

        final Set<ExpectedVariant> expected = Set.of(
                ExpectedVariant.of("X", 34059778, "CATTATAGTAT[X:34062705[", "X", 34062705, "]X:34059778]ATTATAGTATG")
        );
        assertVariants(results, expected);

        final VariantCall call = results.iterator().next();
        assertTrue(call.variantAssemblies()).hasSize(1);
        final AlignedAssembly assembly = call.variantAssemblies().iterator().next().Assembly;
        assertTrue(assembly.Assembly).doesNotContain("TTCCCCCCCCCGCC");
        assertTrue(assembly.Assembly).contains("TTCCCCCGCC");
    }

    private static void assertVariants(final Set<VariantCall> actualRaw, final Set<ExpectedVariant> expected)
    {
        final Set<VariantCall> actual = new HashSet<>(actualRaw); // Make a copy to make debugging easier
        final Map<String, List<VariantCall>> variantsByLeftDescriptor = new HashMap<>();
        for(VariantCall call : actual)
        {
            // Ensure the calls are well-formed, hard-fail at this point.
            assertTrue(call.LeftChromosome).isNotNull();
            assertTrue(call.LeftPosition).isNotEqualTo(0);
            assertTrue(call.LeftDescriptor).isNotNull();

            if(call.RightChromosome == null)
            {
                assertTrue(call.RightPosition).isEqualTo(0);
                assertTrue(call.RightDescriptor).isNull();
            }
            else
            {
                assertTrue(call.RightPosition).isNotEqualTo(0);
                assertTrue(call.RightDescriptor).isNotNull();
            }

            variantsByLeftDescriptor.computeIfAbsent(call.LeftDescriptor, __ -> new ArrayList<>()).add(call);
        }

        final List<String> errors = new ArrayList<>();
        if(actual.size() != expected.size())
            errors.add(String.format("Expected %s variants, but found %s", expected.size(), actual.size()));

        // Did we miss any variants?
        for(ExpectedVariant expectedCall : expected)
        {
            final List<VariantCall> candidates = variantsByLeftDescriptor.get(expectedCall.Left.Descriptor);
            if(candidates == null || candidates.isEmpty())
            {
                errors.add(String.format("Failed to find any candidates for expected variant %s", expectedCall));
                break;
            }

            if(candidates.size() == 1)
            {
                final VariantCall actualCall = candidates.get(0);
                if(!Objects.equals(new VariantSide(actualCall.LeftChromosome, actualCall.LeftPosition, actualCall.LeftDescriptor), expectedCall.Left)
                        || !Objects.equals(new VariantSide(actualCall.RightChromosome, actualCall.RightPosition, actualCall.RightDescriptor), expectedCall.Right))
                    errors.add(String.format("Call %s was not found exactly, found instead: %s", expectedCall, actualCall));
                else
                    actual.remove(actualCall);
            }
            else
            {
                @Nullable
                VariantCall found = null;
                for(VariantCall actualCall : candidates)
                {
                    if(Objects.equals(new VariantSide(actualCall.LeftChromosome, actualCall.LeftPosition, actualCall.LeftDescriptor), expectedCall.Left)
                            && Objects.equals(new VariantSide(actualCall.RightChromosome, actualCall.RightPosition, actualCall.RightDescriptor), expectedCall.Right))
                    {
                        found = actualCall;
                        break;
                    }
                }
                if(found == null)
                    errors.add(String.format("Failed to find match for %s, candidates were: %s", expectedCall, candidates));
                else
                    actual.remove(found);
            }
        }

        // Are there variants left-over that we didn't expect?
        for(VariantCall call : actual)
            errors.add(String.format("Unexpected VariantCall: %s", call));

        assertTrue(errors)
                .withFailMessage(() -> "\n" + String.join("\n", errors))
                .isEmpty();
    }

    @SuppressWarnings("unused") // Call manually to generate expected calls
    private static String createExpectedVariants(final Set<VariantCall> calls)
    {
        if(calls.isEmpty())
            return "final Set<ExpectedVariant> expected = Set.of();";

        final StringBuilder sb = new StringBuilder();
        sb.append("final Set<ExpectedVariant> expected = Set.of(\n");
        for(VariantCall call : calls)
        {
            if(call.isSingleSided())
            {
                sb.append("\t\tExpectedVariant.of(");
                sb.append("\"").append(call.LeftChromosome).append("\", ").append(call.LeftPosition)
                        .append(", \"").append(call.LeftDescriptor).append("\"),\n");
            }
            else
            {
                sb.append("\t\tExpectedVariant.of(\"").append(call.LeftChromosome).append("\", ").append(call.LeftPosition)
                        .append(", \"").append(call.LeftDescriptor).append("\", \"")
                        .append(call.RightChromosome).append("\", ").append(call.RightPosition)
                        .append(", \"").append(call.RightDescriptor).append("\"),\n");
            }
        }
        sb.setLength(sb.length() - 2); // Remove last trailing comma
        sb.append("\n);");

        return sb.toString().replace("\t", "    ");
    }

    private static class ExpectedVariant
    {
        public final VariantSide Left;
        @Nullable
        public final VariantSide Right;

        private ExpectedVariant(final VariantSide left)
        {
            this(left, null);
        }

        private ExpectedVariant(final VariantSide left, @Nullable final VariantSide right)
        {
            Left = left;
            Right = right;
        }

        @Override
        public String toString()
        {
            return Right == null
                    ? String.format("{ %s (SGL) }", Left)
                    : String.format("{ %s -> %s }", Left, Right);
        }

        public static ExpectedVariant of(final String leftChromosome, final int leftPosition, final String leftDescriptor,
                final String rightChromosome, final int rightPosition, final String rightDescriptor)
        {
            return new ExpectedVariant(
                    new VariantSide(leftChromosome, leftPosition, leftDescriptor),
                    new VariantSide(rightChromosome, rightPosition, rightDescriptor)
            );
        }

        public static ExpectedVariant of(final String leftChromosome, final int leftPosition, final String leftDescriptor)
        {
            return new ExpectedVariant(
                    new VariantSide(leftChromosome, leftPosition, leftDescriptor)
            );
        }
    }

    private static class VariantSide
    {
        public final String Chromosome;
        public final int Position;
        public final String Descriptor;

        private VariantSide(final String chromosome, final int position, final String descriptor)
        {
            Chromosome = chromosome;
            Position = position;
            Descriptor = descriptor;
        }

        @Override
        public boolean equals(final Object o)
        {
            if(this == o)
                return true;
            else if(o == null)
                return Chromosome == null && Descriptor == null && Position == 0;
            else if(getClass() != o.getClass())
                return false;

            final VariantSide that = (VariantSide) o;
            return Position == that.Position
                    && Objects.equals(Chromosome, that.Chromosome)
                    && Objects.equals(Descriptor, that.Descriptor);
        }

        @Override
        public int hashCode()
        {
            int result = Chromosome != null ? Chromosome.hashCode() : 0;
            result = 31 * result + Position;
            result = 31 * result + (Descriptor != null ? Descriptor.hashCode() : 0);
            return result;
        }

        @Override
        public String toString()
        {
            return String.format("%s@%s (%s)", Chromosome, Position, Descriptor);
        }
    }

     */
}
