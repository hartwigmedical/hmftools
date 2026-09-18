package com.hartwig.hmftools.virusdetect;

import static java.util.Comparator.comparing;
import static java.util.Comparator.comparingInt;
import static java.util.stream.Collectors.toMap;

import static com.hartwig.hmftools.virusdetect.VirusConstants.REPORTED_MARGINS;

import java.util.ArrayList;
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
                file, ContigStatsColumn.values(), rows, (contigRow, row) ->
                {
                    ContigSelectionResult result = contigRow.contig();
                    ContigStats stat = result.stats();
                    ViralContig contig = result.contig();

                    row.set(ContigStatsColumn.contig, contig.name());
                    row.set(ContigStatsColumn.virus_name, contig.virusName());
                    row.set(ContigStatsColumn.oncology_group, contig.oncologyGroup());
                    row.set(ContigStatsColumn.contig_length, contig.length());
                    row.set(ContigStatsColumn.read_count, stat.readCount());
                    row.set(ContigStatsColumn.multi_align_reads, stat.multiAlignReads());
                    row.set(ContigStatsColumn.origin_clipped_reads, stat.originClippedReads());
                    row.set(ContigStatsColumn.coverage_fraction, stat.coverageFraction());
                    row.set(ContigStatsColumn.read_votes, stat.readVotes());

                    SummaryStats alignPerRead = stat.alignPerRead();
                    row.set(ContigStatsColumn.align_per_read_mean, alignPerRead.mean());
                    row.set(ContigStatsColumn.align_per_read_min, alignPerRead.min());
                    row.set(ContigStatsColumn.align_per_read_p5, alignPerRead.p5());
                    row.set(ContigStatsColumn.align_per_read_p25, alignPerRead.p25());
                    row.set(ContigStatsColumn.align_per_read_p50, alignPerRead.p50());
                    row.set(ContigStatsColumn.align_per_read_p75, alignPerRead.p75());
                    row.set(ContigStatsColumn.align_per_read_p95, alignPerRead.p95());
                    row.set(ContigStatsColumn.align_per_read_max, alignPerRead.max());

                    SummaryStats depth = stat.depth();
                    row.set(ContigStatsColumn.depth_mean, depth.mean());
                    row.set(ContigStatsColumn.depth_min, depth.min());
                    row.set(ContigStatsColumn.depth_p5, depth.p5());
                    row.set(ContigStatsColumn.depth_p25, depth.p25());
                    row.set(ContigStatsColumn.depth_p50, depth.p50());
                    row.set(ContigStatsColumn.depth_p75, depth.p75());
                    row.set(ContigStatsColumn.depth_p95, depth.p95());
                    row.set(ContigStatsColumn.depth_max, depth.max());

                    SummaryStats alignerScore = stat.alignerScore();
                    row.set(ContigStatsColumn.aligner_score_mean, alignerScore.mean());
                    row.set(ContigStatsColumn.aligner_score_min, alignerScore.min());
                    row.set(ContigStatsColumn.aligner_score_p5, alignerScore.p5());
                    row.set(ContigStatsColumn.aligner_score_p25, alignerScore.p25());
                    row.set(ContigStatsColumn.aligner_score_p50, alignerScore.p50());
                    row.set(ContigStatsColumn.aligner_score_p75, alignerScore.p75());
                    row.set(ContigStatsColumn.aligner_score_p95, alignerScore.p95());
                    row.set(ContigStatsColumn.aligner_score_max, alignerScore.max());

                    row.set(ContigStatsColumn.oncology_group_resolution, contigRow.resolution().name());
                    row.set(ContigStatsColumn.oncology_group_outcome, contigRow.outcome().name());
                    row.set(ContigStatsColumn.filter_status, result.filterStatus().name());
                    row.setOrNull(ContigStatsColumn.vote_share_pre_filter, contigRow.preFilterVoteShare());
                    row.setOrNull(ContigStatsColumn.vote_share_post_filter, contigRow.postFilterVoteShare());
                    row.setOrNull(ContigStatsColumn.vote_share_ratio, contigRow.voteShareOfTop());

                    CandidateSelectionResult candidate = result.candidate();
                    row.setOrNull(ContigStatsColumn.votes_rank, candidate == null ? null : String.valueOf(candidate.votesRank()));
                    row.setOrNull(ContigStatsColumn.comparable, candidate == null ? null : String.valueOf(candidate.comparable()));
                    row.setOrNull(ContigStatsColumn.role, candidate == null ? null : candidate.role().name());
                    row.setOrNull(ContigStatsColumn.challenges_ranks, candidate == null ? null : ranks(candidate.challengesRanks()));
                    row.setOrNull(ContigStatsColumn.challenged_by_ranks, candidate == null ? null : ranks(candidate.challengedByRanks()));
                });

        LOGGER.info("wrote {} contig stats to {}", rows.size(), file);
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
                        marginColumns().stream())
                .toList();

        DelimFileWriter.write(
                file, columns, pairs, (pair, row) ->
                {
                    row.set(PairwiseMarginsColumn.oncology_group, pair.subject().oncologyGroup());
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

    private static double votes(List<ContigSelectionResult> contigs)
    {
        return contigs.stream().mapToDouble(contig -> contig.stats().readVotes()).sum();
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

    private static List<String> marginColumns()
    {
        List<String> columns = new ArrayList<>();
        for(int margin : REPORTED_MARGINS)
        {
            columns.add(marginColumn(margin));
        }
        return columns;
    }

    private static String marginColumn(int margin)
    {
        return "reads_m" + margin;
    }

    private static String asString(Object value)
    {
        return value == null ? null : value.toString();
    }

    private record ContigStatsRow(
            String oncologyGroup,
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
        origin_clipped_reads,
        depth_mean,
        depth_min,
        depth_p5,
        depth_p25,
        depth_p50,
        depth_p75,
        depth_p95,
        depth_max,
        align_per_read_mean,
        align_per_read_min,
        align_per_read_p5,
        align_per_read_p25,
        align_per_read_p50,
        align_per_read_p75,
        align_per_read_p95,
        align_per_read_max,
        aligner_score_mean,
        aligner_score_min,
        aligner_score_p5,
        aligner_score_p25,
        aligner_score_p50,
        aligner_score_p75,
        aligner_score_p95,
        aligner_score_max
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
}
