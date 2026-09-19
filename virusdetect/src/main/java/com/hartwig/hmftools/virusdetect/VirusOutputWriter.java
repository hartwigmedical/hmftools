package com.hartwig.hmftools.virusdetect;

import static java.util.Comparator.comparing;
import static java.util.Comparator.comparingInt;
import static java.util.stream.Collectors.toMap;

import static com.hartwig.hmftools.virusdetect.VirusConstants.REPORTED_MARGINS;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

public class VirusOutputWriter
{
    private static final Logger LOGGER = LogManager.getLogger(VirusOutputWriter.class);

    // One row per contig with any alignment: its statistics, and how representative selection judged it.
    public static void writeContigStats(String file, List<OncologyGroupSelection> selections)
    {
        List<ContigStatsRow> rows = selections.stream()
                .flatMap(VirusOutputWriter::contigStatsRows)
                .sorted(comparingInt((ContigStatsRow row) -> row.contig().stats().readCount()).reversed()
                        .thenComparing(row -> row.contig().contig().name()))
                .toList();

        DelimFileWriter.write(
                file, CONTIG_STATS_COLUMNS, rows, (contigRow, row) ->
                {
                    ContigSelectionResult result = contigRow.contig();
                    ContigSupport stat = result.stats();
                    ViralContig contig = result.contig();

                    row.set(ContigStatsColumn.contig, contig.name());
                    row.set(ContigStatsColumn.virus_name, contig.virusName());
                    row.set(ContigStatsColumn.oncology_group, contig.oncologyGroup().name());
                    row.set(ContigStatsColumn.contig_length, contig.length());
                    row.set(ContigStatsColumn.read_count, stat.readCount());
                    row.set(ContigStatsColumn.multi_align_reads, stat.multiAlignReads());
                    row.set(ContigStatsColumn.origin_clipped_reads, stat.originClippedReads());
                    row.set(ContigStatsColumn.coverage_fraction, stat.coverageFraction());
                    row.set(ContigStatsColumn.read_votes, stat.readVotes());

                    writeSummaryStats(row, DEPTH_STATS_COLUMNS, stat.depth());
                    writeSummaryStats(row, ALIGN_PER_READ_STATS_COLUMNS, stat.alignPerRead());
                    writeSummaryStats(row, ALIGNER_SCORE_STATS_COLUMNS, stat.alignerScore());

                    row.set(ContigStatsColumn.oncology_group_resolution, contigRow.resolution().name());
                    row.set(ContigStatsColumn.oncology_group_outcome, contigRow.outcome().name());
                    row.set(ContigStatsColumn.filter_status, result.filterStatus().name());
                    row.setOrNull(ContigStatsColumn.vote_share_pre_filter, contigRow.preFilterVoteShare());
                    row.setOrNull(ContigStatsColumn.vote_share_post_filter, contigRow.postFilterVoteShare());
                    row.setOrNull(ContigStatsColumn.vote_share_ratio, contigRow.voteShareOfTop());

                    // Note unset columns are written as null.
                    CandidateSelectionResult candidate = result.candidate();
                    if(candidate != null)
                    {
                        row.set(ContigStatsColumn.votes_rank, candidate.votesRank());
                        row.set(ContigStatsColumn.comparable, candidate.comparable());
                        row.set(ContigStatsColumn.role, candidate.role().name());
                        row.set(ContigStatsColumn.challenges_ranks, ranks(candidate.challengesRanks()));
                        row.set(ContigStatsColumn.challenged_by_ranks, ranks(candidate.challengedByRanks()));
                    }
                });

        LOGGER.info("wrote {} contig stats to {}", rows.size(), file);
    }

    private static final String DEPTH_STATS_COLUMNS = "depth";
    private static final String ALIGN_PER_READ_STATS_COLUMNS = "align_per_read";
    private static final String ALIGNER_SCORE_STATS_COLUMNS = "aligner_score";
    private static final List<String> CONTIG_STATS_COLUMNS =
            Stream.concat(
                            Stream.of(ContigStatsColumn.values()).map(Enum::name),
                            Stream.of(DEPTH_STATS_COLUMNS, ALIGN_PER_READ_STATS_COLUMNS, ALIGNER_SCORE_STATS_COLUMNS)
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
        return group + "_" + field;
    }

    private static double votes(List<ContigSelectionResult> contigs)
    {
        return contigs.stream().mapToDouble(contig -> contig.stats().readVotes()).sum();
    }

    private record ContigStatsRow(
            OncologyGroup oncologyGroup,
            OncologyGroupOutcome outcome,
            ContigSelectionResult contig,
            @Nullable Double preFilterVoteShare,
            @Nullable Double postFilterVoteShare,
            @Nullable Double voteShareOfTop
    )
    {
        private OncologyGroupResolution resolution()
        {
            return outcome.resolution();
        }
    }

    private static Stream<ContigStatsRow> contigStatsRows(OncologyGroupSelection selection)
    {
        double allContigVotes = votes(selection.contigs());
        double candidateVotes = votes(selection.candidates());
        double topCandidateVotes = selection.candidates().stream()
                .mapToDouble(candidate -> candidate.stats().readVotes()).max().orElse(0.0);

        return selection.contigs().stream().map(contig ->
        {
            double contigVotes = contig.stats().readVotes();
            return new ContigStatsRow(
                    selection.oncologyGroup(), selection.outcome(), contig,
                    share(contigVotes, allContigVotes), share(contigVotes, candidateVotes),
                    share(contigVotes, topCandidateVotes));
        });
    }

    private enum ContigStatsColumn
    {
        contig,
        virus_name,
        oncology_group,
        oncology_group_resolution,
        oncology_group_outcome,
        role,
        filter_status,
        contig_length,
        coverage_fraction,
        read_count,
        read_votes,
        votes_rank,
        vote_share_pre_filter,
        vote_share_post_filter,
        vote_share_ratio,
        comparable,
        challenges_ranks,
        challenged_by_ranks,
        multi_align_reads,
        origin_clipped_reads
    }

    @Nullable
    private static Double share(double votes, double total)
    {
        return total > 0 ? votes / total : null;
    }

    private static String ranks(List<Integer> ranks)
    {
        return ranks.stream().map(String::valueOf).collect(Collectors.joining(","));
    }

    // One row per ordered within-oncology-group contig pair: the reads they share, and how many of those fit the subject
    // better by at least each reported margin. Includes prefiltered contigs. Verbose/debug only, for tuning.
    public static void writePairwiseMargins(String file, PairwiseMargins margins, List<OncologyGroupSelection> selections)
    {
        Map<ViralContig, Integer> votesRankByContig = selections.stream()
                .flatMap(selection -> selection.candidates().stream())
                .collect(toMap(ContigSelectionResult::contig, result -> result.candidate().votesRank()));

        List<PairwiseMargins.ContigPair> pairs = margins.pairs().stream()
                .sorted(comparing((PairwiseMargins.ContigPair pair) -> pair.subject().name())
                        .thenComparing(pair -> pair.opponent().name()))
                .toList();

        List<String> columns = Stream.concat(
                        Stream.of(PairwiseMarginsColumn.values()).map(Enum::name),
                        MARGIN_COLUMNS.stream())
                .toList();

        DelimFileWriter.write(
                file, columns, pairs, (pair, row) ->
                {
                    row.set(PairwiseMarginsColumn.oncology_group, pair.subject().oncologyGroup().name());
                    row.set(PairwiseMarginsColumn.subject_contig, pair.subject().name());
                    row.setOrNull(PairwiseMarginsColumn.subject_rank, asString(votesRankByContig.get(pair.subject())));
                    row.set(PairwiseMarginsColumn.opponent_contig, pair.opponent().name());
                    row.setOrNull(PairwiseMarginsColumn.opponent_rank, asString(votesRankByContig.get(pair.opponent())));
                    row.set(PairwiseMarginsColumn.shared_reads, margins.sharedReads(pair.subject(), pair.opponent()));
                    for(int margin : REPORTED_MARGINS)
                    {
                        row.set(marginColumn(margin), margins.readsWinningBy(pair.subject(), pair.opponent(), margin));
                    }
                });

        LOGGER.info("wrote {} pairwise margin rows to {}", pairs.size(), file);
    }

    private enum PairwiseMarginsColumn
    {
        oncology_group,
        subject_contig,
        subject_rank,
        opponent_contig,
        opponent_rank,
        shared_reads
    }

    private static final List<String> MARGIN_COLUMNS = REPORTED_MARGINS.stream().map(VirusOutputWriter::marginColumn).toList();

    private static String marginColumn(int margin)
    {
        return "reads_m" + margin;
    }

    private static String asString(Object value)
    {
        return value == null ? null : value.toString();
    }
}
