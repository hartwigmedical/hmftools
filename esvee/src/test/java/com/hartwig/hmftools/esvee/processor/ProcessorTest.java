package com.hartwig.hmftools.esvee.processor;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.Junction;

import org.immutables.value.Value;
import org.jetbrains.annotations.Nullable;
import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class ProcessorTest
{
    public static final String FULL_SOMATIC_BAM = "/Users/james/code/data/COLO829/COLO829v003T.bam";
    public static final String FULL_GERMLINE_BAM = "/Users/james/code/data/COLO829/COLO829v003R.bam";

    public static final String SOMATIC_BAM = "/Users/james/code/data/COLO829/sv-prep2/COLO829v003T.sv_prep_sorted.bam";
    public static final String GERMLINE_BAM = "/Users/james/code/data/COLO829/sv-prep2/COLO829v003R.sv_prep_sorted.bam";
    public static final String JUNCTION_FILE = "/Users/james/code/data/COLO829/sv-prep2/COLO829v003T.sv_prep.junctions.tsv";

    private static final String REFERENCE_GENOME = "/Users/james/code/data/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta";
    private static final String REFERENCE_GENOME_INDEX = "/Users/james/code/data/refgenome/Homo_sapiens.GRCh37.GATK.illumina.img";

    private int containsNearJunction(final Map<String, Set<Integer>> positions, final Junction junction, final int tolerance)
    {
        return containsNearJunction(positions, junction.chromosome(), junction.position(), tolerance);
    }

    private int containsNearJunction(final Map<String, Set<Integer>> positions, final String chromosome, final int position, final int tolerance)
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


    @Value.Immutable
    public interface Expected
    {
        String curation();
        String chromosome();
        int position();
        String assemblies();
    }

    /*

    @Test
    public void processTruePositives() throws IOException
    {
        final SVAConfig config = HMFConfig.load(Map.of(
                        //"bam_file", "/Users/james/code/data/run/old/tumor.sv_prep.sorted.bam",
                        //"germline_bam", "/Users/james/code/data/run/old/ref.sv_prep.sorted.bam",
                        "bam_file", SOMATIC_BAM,
                        "germline_bam", GERMLINE_BAM,
                        "output_file", "output.vcf",
                        "ref_genome", "/Users/james/code/data/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta",
                        "ref_genome_index", "/Users/james/code/data/refgenome/Homo_sapiens.GRCh37.GATK.illumina.img",
                        "junction_file", JUNCTION_FILE,
                        "create_html_summaries", "true"
                        , "output_bam_file", "output.bam"
                        //,"debug", "true"
                        //,"timeouts_enabled", "false"
                                , "threads", "1"
                ),
                SVAConfig.class, ImmutableSVAConfig.builder());
        final Context context = Context.create(config);
        final Processor processor = new Processor(context);

        final List<Expected> expected = new CSVReader<Expected>(ImmutableExpected.class, "expected.csv").readToEnd();
        final Map<String, Set<Integer>> toProcessPositions = new HashMap<>();
        for (final Expected e : expected)
            toProcessPositions.computeIfAbsent(e.chromosome(), ignored -> new HashSet<>()).add(e.position());

        final List<Junction> junctions = JunctionReader.readJunctionFile(new File("/Users/james/code/data/run/old/tumor.sv_prep.junctions.tsv"))
                .filter(junction -> containsNearJunction(toProcessPositions, junction, 5) != -1)
                .collect(Collectors.toList());

        final List<VariantCall> results = processor.run(junctions);

        final Map<String, Map<Integer, List<VariantCall>>> resultsByPosition = new HashMap<>();
        for (final VariantCall call : results)
        {
            if (call.LeftChromosome != null)
                resultsByPosition.computeIfAbsent(call.LeftChromosome, ignored -> new HashMap<>())
                        .computeIfAbsent(containsNearJunction(toProcessPositions, call.LeftChromosome, call.LeftPosition, 5), k -> new ArrayList<>())
                        .add(call);
            if (call.RightChromosome != null)
                resultsByPosition.computeIfAbsent(call.RightChromosome, ignored -> new HashMap<>())
                        .computeIfAbsent(containsNearJunction(toProcessPositions, call.RightChromosome, call.RightPosition, 5), k -> new ArrayList<>())
                        .add(call);
        }

        final List<ImmutableExpected> updatedExpected = new ArrayList<>();
        for (final Expected e : expected)
        {
            final List<VariantCall> resultsForPosition = resultsByPosition.getOrDefault(e.chromosome(), Map.of())
                    .getOrDefault(e.position(), List.of());
            final String assemblyCounts = resultsForPosition.stream().map(r -> String.valueOf(r.associatedAssemblies().size())).collect(Collectors.joining(","));
            updatedExpected.add(ImmutableExpected.builder()
                    .curation(e.curation())
                    .chromosome(e.chromosome())
                    .position(e.position())
                    .assemblies(assemblyCounts)
                    .build());
        }
        updatedExpected.sort(Comparator.comparing(ImmutableExpected::chromosome)
                .thenComparing(ImmutableExpected::position));
        CSVWriter.writeCSV("actual.csv", ImmutableExpected.class, updatedExpected);
    }

    @Test
    public void processSingle()
    {
        final List<Pair<String, Integer>> locations = List.of(
                //Pair.of("9", 28031835),
                //Pair.of("8", 131248116),
                //Pair.of("1", 5447171) // Poor mappability
                //Pair.of("3", 145410927)
                //Pair.of("19", 54509005)
                //Pair.of("1", 210261697),
                //Pair.of("1", 210262333)
                //Pair.of("1", 821283)
                //Pair.of("11", 76903948), // Timeout during primary assembly
                //Pair.of("9", 30338965) // Phantom extension
                //Pair.of("1", 81787809) // Missed feature?
                //Pair.of("3", 184701983) // Timeout, legit insert

                //                Pair.of("3", 26663922),
                //                Pair.of("3", 26664498)

                //                Pair.of("2", 113096430),
                //                Pair.of("2", 113096447),
                //                Pair.of("2", 113096449),
                //                Pair.of("2", 113096473)

                //                Pair.of("9", 28031863),
                //                Pair.of("9", 28034467),
                //                Pair.of("9", 28031835),
                //                Pair.of("9", 28059140)

                //                Pair.of("4", 79170203),
                //                Pair.of("4", 79170205),
                //                Pair.of("4", 79170207),
                //                Pair.of("4", 79170217),
                //                Pair.of("4", 79170219),
                //                Pair.of("4", 79170234)

                //Pair.of("1", 4204544)
                //Pair.of("11", 39948209),
                //Pair.of("11", 39948217)

                //Pair.of("9", 30338965) // TODO:

                //
                //                Pair.of("11", 720944),
                //                Pair.of("11", 721008)

                //                Pair.of("9", 30532655),
                //                Pair.of("9", 30532660),
                //                Pair.of("9", 30532665),
                //                Pair.of("9", 30532673)
        );
        final SVAConfig config = HMFConfig.load(Map.of(
                        "bam_file", SOMATIC_BAM,
                        "germline_bam", GERMLINE_BAM,
                        //                        "bam_file", FULL_SOMATIC_BAM,
                        //                        "germline_bam", FULL_GERMLINE_BAM,
                        "output_file", "output.vcf",
                        "ref_genome", REFERENCE_GENOME,
                        "ref_genome_index", "/Users/james/code/data/refgenome/Homo_sapiens.GRCh37.GATK.illumina.img",
                        "junction_file", JUNCTION_FILE,
                        "timeouts_enabled", "false",
                        "debug", "true"
                        //"drop_germline", "true"
                        , "create_html_summaries", "true"
                        , "threads", "1"
                ),
                SVAConfig.class, ImmutableSVAConfig.builder());
        final Context context = Context.create(config);
        final Processor processor = new Processor(context);

        final Map<String, Set<Integer>> toProcessPositions = new HashMap<>();
        for (final var pair : locations)
            toProcessPositions.computeIfAbsent(pair.getKey(), ignored -> new HashSet<>()).add(pair.getValue());

        final List<Junction> junctions = JunctionReader.readJunctionFile(config.junctionFile())
                .filter(junction -> containsNearJunction(toProcessPositions, junction, 5) != -1)
                .collect(Collectors.toList());

        processor.run(junctions);
    }

    @Test
    public void processProblems()
    {
        final List<Pair<String, Integer>> locations = List.of(
                Pair.of("20", 62728305),
                Pair.of("20", 62728393),
                Pair.of("20", 62728406),
                Pair.of("20", 62728394),
                Pair.of("20", 62728389)
        );
        final SVAConfig config = HMFConfig.load(Map.of(
                        "bam_file", SOMATIC_BAM,
                        "germline_bam", GERMLINE_BAM,
                        "output_file", "output.vcf",
                        "ref_genome", REFERENCE_GENOME,
                        "ref_genome_index", REFERENCE_GENOME_INDEX,
                        "junction_file", JUNCTION_FILE,
                        "timeouts_enabled", "false",
                        "debug", "true"
                        , "create_html_summaries", "true"
                ),
                SVAConfig.class, ImmutableSVAConfig.builder());
        final Context context = Context.create(config);
        final Processor processor = new Processor(context);

        final Map<String, Set<Integer>> toProcessPositions = new HashMap<>();
        for (final var pair : locations)
            toProcessPositions.computeIfAbsent(pair.getKey(), ignored -> new HashSet<>()).add(pair.getValue());

        final List<Junction> junctions = JunctionReader.readJunctionFile(config.junctionFile())
                .filter(junction -> containsNearJunction(toProcessPositions, junction, 5) != -1)
                .collect(Collectors.toList());

        processor.run(junctions);
    }

    @Test
    public void processAll()
    {
        final SVAConfig config = HMFConfig.load(Map.of(
                        "bam_file", SOMATIC_BAM,
                        "germline_bam", GERMLINE_BAM,
                        "output_file", "output.vcf",
                        "ref_genome", "/Users/james/code/data/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta",
                        "ref_genome_index", "/Users/james/code/data/refgenome/Homo_sapiens.GRCh37.GATK.illumina.img",
                        "junction_file", JUNCTION_FILE,
                        "create_html_summaries", "true",
                        "create_diagrams", "false"
                        , "drop_germline", "true"
                        , "html_summaries_folder", "SummariesFull"
                ),
                SVAConfig.class, ImmutableSVAConfig.builder());
        assert config != null;
        final Context context = Context.create(config);
        final Processor processor = new Processor(context);

        processor.run();
    }
     */
}