package com.hartwig.hmftools.viridian.app;

import static java.util.Comparator.comparing;
import static java.util.Comparator.reverseOrder;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.joining;
import static java.util.stream.Collectors.summingDouble;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.viridian.integration.Integration;
import com.hartwig.hmftools.viridian.integration.align.ViralInsertAlignment;
import com.hartwig.hmftools.viridian.integration.extract.BreakendSupport;
import com.hartwig.hmftools.viridian.integration.extract.CandidateSv;
import com.hartwig.hmftools.viridian.integration.extract.HostBreakend;
import com.hartwig.hmftools.viridian.integration.extract.InsertRepeat;
import com.hartwig.hmftools.viridian.profiling.VirusProfile;
import com.hartwig.hmftools.viridian.profiling.common.GenomeStats;
import com.hartwig.hmftools.viridian.profiling.common.SummaryStats;
import com.hartwig.hmftools.viridian.profiling.support.GenomeSupport;
import com.hartwig.hmftools.viridian.reference.VirusGenome;
import com.hartwig.hmftools.viridian.reference.VirusGroup;

import org.jetbrains.annotations.Nullable;

public class OutputWriter
{
    public static void writeGenomeSupport(
            String file, List<GenomeSupport> genomeSupports, List<VirusGenome> representatives)
    {
        Set<VirusGenome> representativeGenomes = Set.copyOf(representatives);
        Map<VirusGroup, Double> groupVotes = genomeSupports.stream()
                .collect(groupingBy(support -> support.genome().virusGroup(), summingDouble(GenomeSupport::readVotes)));
        List<GenomeSupport> rows = genomeSupports.stream()
                .sorted(comparing((GenomeSupport support) -> support.genome().virusGroup())
                        .thenComparing(GenomeSupport::readVotes, reverseOrder())
                        .thenComparing(GenomeSupport::genome))
                .toList();

        DelimFileWriter.write(
                file, GENOME_SUPPORT_COLUMNS, rows, (support, row) ->
                {
                    GenomeStats stats = support.stats();
                    VirusGenome genome = support.genome();

                    row.set(GenomeSupportColumn.contig, genome.contigName());
                    row.set(GenomeSupportColumn.virusName, genome.virus().name());
                    row.set(GenomeSupportColumn.virusGroup, genome.virusGroup().name());
                    row.set(GenomeSupportColumn.genomeLength, genome.length());
                    row.set(GenomeSupportColumn.alignedReads, stats.readCount());
                    row.set(GenomeSupportColumn.multiAlignReads, support.multiAlignReads());
                    row.set(GenomeSupportColumn.originClippedReads, stats.originClippedReads());
                    row.set(GenomeSupportColumn.coverage, stats.coverageFraction());
                    row.set(GenomeSupportColumn.readVotes, support.readVotes());

                    writeSummaryStats(row, DEPTH_STATS_COLUMNS, stats.depth());
                    writeSummaryStats(row, ALIGN_PER_READ_STATS_COLUMNS, support.alignmentsPerRead());
                    writeSummaryStats(row, ALIGNER_SCORE_STATS_COLUMNS, stats.alignerScore());

                    row.setOrNull(
                            GenomeSupportColumn.voteFraction,
                            fraction(support.readVotes(), groupVotes.get(genome.virusGroup())));
                    row.set(GenomeSupportColumn.topMatchReads, support.topMatchReads());
                    row.set(GenomeSupportColumn.isRepresentative, representativeGenomes.contains(genome));
                });
    }

    private enum GenomeSupportColumn
    {
        contig,
        virusName,
        virusGroup,
        isRepresentative,
        genomeLength,
        coverage,
        alignedReads,
        readVotes,
        voteFraction,
        topMatchReads,
        multiAlignReads,
        originClippedReads
    }

    private static final String DEPTH_STATS_COLUMNS = "depth";
    private static final String ALIGNER_SCORE_STATS_COLUMNS = "alignerScore";

    private static final String ALIGN_PER_READ_STATS_COLUMNS = "alignPerRead";

    private static final List<String> GENOME_SUPPORT_COLUMNS =
            Stream.concat(
                            Stream.of(GenomeSupportColumn.values()).map(Enum::name),
                            Stream.of(DEPTH_STATS_COLUMNS, ALIGN_PER_READ_STATS_COLUMNS, ALIGNER_SCORE_STATS_COLUMNS)
                                    .flatMap(group -> SummaryStats.FIELD_NAMES.stream()
                                            .map(field -> summaryStatsColumn(group, field))))
                    .toList();

    public static void writeVirusProfiles(String file, List<VirusProfile> profiles)
    {
        List<VirusProfile> rows = profiles.stream().sorted(comparing(VirusProfile::virusGroup)).toList();

        DelimFileWriter.write(
                file, VIRUS_PROFILE_COLUMNS, rows, (profile, row) ->
                {
                    row.set(VirusProfileColumn.virusGroup, profile.virusGroup().name());

                    GenomeStats stats = profile.representativeStats();
                    VirusGenome genome = stats.genome();
                    row.set(VirusProfileColumn.contig, genome.contigName());
                    row.set(VirusProfileColumn.virusName, genome.virus().name());
                    row.set(VirusProfileColumn.genomeLength, genome.length());
                    row.set(VirusProfileColumn.alignedReads, stats.readCount());
                    row.set(VirusProfileColumn.originClippedReads, stats.originClippedReads());
                    row.set(VirusProfileColumn.coverage, stats.coverageFraction());

                    List<VirusGenome> secondaries = profile.secondaryGenomes();
                    row.setOrNull(VirusProfileColumn.secondaryContigs, joinOrNull(secondaries, VirusGenome::contigName));
                    row.setOrNull(
                            VirusProfileColumn.secondaryVirusNames,
                            joinOrNull(secondaries, secondary -> secondary.virus().name()));

                    writeSummaryStats(row, DEPTH_STATS_COLUMNS, stats.depth());
                    writeSummaryStats(row, ALIGNER_SCORE_STATS_COLUMNS, stats.alignerScore());
                });
    }

    private enum VirusProfileColumn
    {
        virusGroup,
        virusName,
        contig,
        genomeLength,
        coverage,
        alignedReads,
        originClippedReads,
        secondaryContigs,
        secondaryVirusNames,
    }

    private static final List<String> VIRUS_PROFILE_COLUMNS =
            Stream.concat(
                            Stream.of(VirusProfileColumn.values()).map(Enum::name),
                            Stream.of(DEPTH_STATS_COLUMNS, ALIGNER_SCORE_STATS_COLUMNS)
                                    .flatMap(group -> SummaryStats.FIELD_NAMES.stream()
                                            .map(field -> summaryStatsColumn(group, field))))
                    .toList();

    private static void writeSummaryStats(DelimFileWriter.Row row, String group, @Nullable SummaryStats stats)
    {
        if(stats == null)
        {
            return;
        }

        stats.fieldValues().forEach((field, value) -> row.set(summaryStatsColumn(group, field), value));
    }

    private static String summaryStatsColumn(String group, String field)
    {
        return group + Character.toUpperCase(field.charAt(0)) + field.substring(1);
    }

    @Nullable
    private static Double fraction(double votes, double total)
    {
        return total > 0 ? votes / total : null;
    }

    @Nullable
    private static String joinOrNull(List<VirusGenome> genomes, Function<VirusGenome, String> field)
    {
        return genomes.isEmpty() ? null : genomes.stream().map(field).collect(joining(ITEM_DELIM));
    }

    public static void writeIntegrations(String file, List<Integration> integrations)
    {
        DelimFileWriter.write(
                file, IntegrationColumn.values(), integrations, (integration, row) ->
                {
                    CandidateSv sv = integration.sv();
                    row.set(IntegrationColumn.svType, sv.type().name());
                    row.set(IntegrationColumn.filter, sv.filter());

                    HostBreakend startBreakend = sv.startBreakend();
                    BreakendSupport startSupport = startBreakend.support();
                    row.set(IntegrationColumn.startId, startBreakend.id());
                    row.set(IntegrationColumn.startChromosome, startBreakend.basePosition().Chromosome);
                    row.set(IntegrationColumn.startPosition, startBreakend.basePosition().Position);
                    row.set(IntegrationColumn.startOrientation, startBreakend.orientation().asByte());
                    row.set(IntegrationColumn.startTumorVariantFrags, startSupport.tumorVariantFragments());
                    row.set(IntegrationColumn.startTumorRefFrags, startSupport.tumorReferenceFragments());
                    row.set(IntegrationColumn.startTumorAf, startSupport.tumorAlleleFrequency());
                    row.setOrNull(
                            IntegrationColumn.startNormalVariantFrags.name(),
                            startSupport.normalVariantFragments());
                    row.setOrNull(
                            IntegrationColumn.startNormalRefFrags.name(),
                            startSupport.normalReferenceFragments());

                    // Note unset columns are written as null, leaving a SGL's second breakend blank.
                    HostBreakend endBreakend = sv.endBreakend();
                    if(endBreakend != null)
                    {
                        BreakendSupport endSupport = endBreakend.support();
                        row.set(IntegrationColumn.endId, endBreakend.id());
                        row.set(IntegrationColumn.endChromosome, endBreakend.basePosition().Chromosome);
                        row.set(IntegrationColumn.endPosition, endBreakend.basePosition().Position);
                        row.set(IntegrationColumn.endOrientation, endBreakend.orientation().asByte());
                        row.set(IntegrationColumn.endTumorVariantFrags, endSupport.tumorVariantFragments());
                        row.set(IntegrationColumn.endTumorRefFrags, endSupport.tumorReferenceFragments());
                        row.set(IntegrationColumn.endTumorAf, endSupport.tumorAlleleFrequency());
                        row.setOrNull(
                                IntegrationColumn.endNormalVariantFrags.name(),
                                endSupport.normalVariantFragments());
                        row.setOrNull(
                                IntegrationColumn.endNormalRefFrags.name(),
                                endSupport.normalReferenceFragments());
                    }

                    row.set(IntegrationColumn.lineSite, sv.lineSite());
                    InsertRepeat insertRepeat = sv.insertRepeat();
                    if(insertRepeat != null)
                    {
                        row.set(IntegrationColumn.repeatClass, insertRepeat.repeatClass());
                        row.set(IntegrationColumn.repeatType, insertRepeat.repeatType());
                        row.set(IntegrationColumn.repeatCoverage, insertRepeat.coverage());
                    }
                    row.set(IntegrationColumn.insertHostAlignments, sv.insertHostAlignments());
                    row.set(IntegrationColumn.insertSeqLength, sv.insertSequence().length());
                    row.set(IntegrationColumn.insertSequence, sv.insertSequence());

                    ViralInsertAlignment alignment = integration.alignment();
                    VirusGenome genome = alignment.genome();
                    row.set(IntegrationColumn.virusName, genome.virus().name());
                    row.set(IntegrationColumn.virusGroup, genome.virusGroup().name());
                    row.set(IntegrationColumn.virusContig, genome.contigName());
                    row.set(IntegrationColumn.virusPosition, alignment.position());
                    row.set(IntegrationColumn.virusOrientation, alignment.orientation().asByte());
                    row.set(IntegrationColumn.alignCigar, alignment.cigar().toString());
                    row.set(IntegrationColumn.alignedLength, alignment.alignedLength());
                    row.set(IntegrationColumn.alignerScore, alignment.alignerScore());
                    row.set(IntegrationColumn.alignedEditDistance, alignment.alignedEditDistance());
                });
    }

    private enum IntegrationColumn
    {
        svType,
        filter,
        startId,
        startChromosome,
        startPosition,
        startOrientation,
        endId,
        endChromosome,
        endPosition,
        endOrientation,
        startTumorVariantFrags,
        startTumorRefFrags,
        startTumorAf,
        startNormalVariantFrags,
        startNormalRefFrags,
        endTumorVariantFrags,
        endTumorRefFrags,
        endTumorAf,
        endNormalVariantFrags,
        endNormalRefFrags,
        lineSite,
        repeatClass,
        repeatType,
        repeatCoverage,
        insertHostAlignments,
        insertSeqLength,
        insertSequence,
        virusName,
        virusGroup,
        virusContig,
        virusPosition,
        virusOrientation,
        alignCigar,
        alignedLength,
        alignerScore,
        alignedEditDistance
    }

    public static void writeRunManifest(String file, ViridianConfig config)
    {
        DelimFileWriter.write(
                file, ManifestColumn.values(), List.of(config), (runConfig, row) ->
                {
                    row.set(ManifestColumn.candidateReadSource, runConfig.candidateReadSource().name());
                    row.set(
                            ManifestColumn.extraVirusDecoyContigs,
                            String.join(ITEM_DELIM, runConfig.extraVirusDecoyContigs().stream().sorted().toList()));
                    row.set(ManifestColumn.kmerFilter, runConfig.kmerFilterEnabled());
                });
    }

    private enum ManifestColumn
    {
        candidateReadSource,
        extraVirusDecoyContigs,
        kmerFilter
    }
}
