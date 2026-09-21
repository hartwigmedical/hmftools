package com.hartwig.hmftools.viridian.app;

import static java.util.Comparator.comparing;
import static java.util.Comparator.comparingInt;
import static java.util.stream.Collectors.toMap;

import static com.hartwig.hmftools.viridian.common.ViridianConstants.REPORTED_MARGINS;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.viridian.common.SummaryStats;
import com.hartwig.hmftools.viridian.detection.contig_support.ContigSupport;
import com.hartwig.hmftools.viridian.reference.ViralContig;
import com.hartwig.hmftools.viridian.selection.OncologyGroupRepresentativeSelection;
import com.hartwig.hmftools.viridian.selection.PairwiseMargins;
import com.hartwig.hmftools.viridian.selection.RepresentativeContigCandidate;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

public class OutputWriter
{
    private static final Logger LOGGER = LogManager.getLogger(OutputWriter.class);

    // TODO: write in order: oncology group (lexicographically, ascending), read votes (descending)
    public static void writeContigInfo(String file, List<OncologyGroupRepresentativeSelection> selections)
    {
        Map<ViralContig, Integer> votesRanks = candidateVotesRanks(selections);
        List<ContigInfoRow> rows = selections.stream()
                .flatMap(OutputWriter::contigInfoRows)
                .sorted(comparingInt((ContigInfoRow row) -> row.support().readCount()).reversed()
                        .thenComparing(row -> row.support().contig().name()))
                .toList();

        DelimFileWriter.write(
                file, CONTIG_INFO_COLUMNS, rows, (contigRow, row) ->
                {
                    ContigSupport support = contigRow.support();
                    ViralContig contig = support.contig();

                    row.set(ContigInfoColumn.contig, contig.name());
                    row.set(ContigInfoColumn.virus_name, contig.virusName());
                    row.set(ContigInfoColumn.oncology_group, contig.oncologyGroup().name());
                    row.set(ContigInfoColumn.contig_length, contig.length());
                    row.set(ContigInfoColumn.read_count, support.readCount());
                    row.set(ContigInfoColumn.multi_align_reads, support.multiAlignReads());
                    row.set(ContigInfoColumn.origin_clipped_reads, support.originClippedReads());
                    row.set(ContigInfoColumn.coverage_fraction, support.coverageFraction());
                    row.set(ContigInfoColumn.read_votes, support.readVotes());

                    writeSummaryStats(row, DEPTH_STATS_COLUMNS, support.depth());
                    writeSummaryStats(row, ALIGN_PER_READ_STATS_COLUMNS, support.alignmentsPerRead());
                    writeSummaryStats(row, ALIGNER_SCORE_STATS_COLUMNS, support.alignerScore());

                    row.set(ContigInfoColumn.oncology_group_resolution, contigRow.selection().resolution().name());
                    row.set(ContigInfoColumn.oncology_group_outcome, contigRow.selection().outcome().name());
                    row.set(ContigInfoColumn.filter_status, support.filterStatus().name());
                    row.setOrNull(ContigInfoColumn.vote_share_pre_filter, contigRow.preFilterVoteShare());
                    row.setOrNull(ContigInfoColumn.vote_share_post_filter, contigRow.postFilterVoteShare());
                    row.setOrNull(ContigInfoColumn.vote_share_ratio, contigRow.voteShareOfTop());

                    // Note unset columns are written as null.
                    RepresentativeContigCandidate candidate = contigRow.candidate();
                    if(candidate != null)
                    {
                        row.set(ContigInfoColumn.votes_rank, candidate.votesRank());
                        row.set(ContigInfoColumn.comparable, candidate.comparable());
                        row.set(ContigInfoColumn.role, candidate.role().name());
                        row.set(ContigInfoColumn.challenges_ranks, formatRanks(votesRanks, candidate.challenges()));
                        row.set(ContigInfoColumn.challenged_by_ranks, formatRanks(votesRanks, candidate.challengedBy()));
                    }
                });

        LOGGER.info("wrote {} contig info rows to {}", rows.size(), file);
    }

    private static final String DEPTH_STATS_COLUMNS = "depth";
    private static final String ALIGN_PER_READ_STATS_COLUMNS = "align_per_read";
    private static final String ALIGNER_SCORE_STATS_COLUMNS = "aligner_score";
    private static final List<String> CONTIG_INFO_COLUMNS =
            Stream.concat(
                            Stream.of(ContigInfoColumn.values()).map(Enum::name),
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

    private record ContigInfoRow(
            OncologyGroupRepresentativeSelection selection,
            ContigSupport support,
            // Null for a contig the support filter rejected, leaving the candidate columns blank.
            @Nullable RepresentativeContigCandidate candidate,
            @Nullable Double preFilterVoteShare,
            @Nullable Double postFilterVoteShare,
            @Nullable Double voteShareOfTop
    )
    {
    }

    private static Stream<ContigInfoRow> contigInfoRows(OncologyGroupRepresentativeSelection selection)
    {
        double candidateVotes = selection.candidates().stream().mapToDouble(c -> c.support().readVotes()).sum();
        double rejectedVotes = selection.rejected().stream().mapToDouble(ContigSupport::readVotes).sum();
        double topCandidateVotes = selection.candidates().stream()
                .mapToDouble(candidate -> candidate.support().readVotes()).max().orElse(0.0);

        Stream<ContigInfoRow> candidates = selection.candidates().stream()
                .map(candidate -> row(
                        selection, candidate.support(), candidate, candidateVotes + rejectedVotes,
                        candidateVotes, topCandidateVotes));
        Stream<ContigInfoRow> rejected = selection.rejected().stream()
                .map(support -> row(
                        selection, support, null, candidateVotes + rejectedVotes,
                        candidateVotes, topCandidateVotes));

        return Stream.concat(candidates, rejected);
    }

    private static ContigInfoRow row(
            OncologyGroupRepresentativeSelection selection, ContigSupport support, @Nullable RepresentativeContigCandidate candidate,
            double allVotes, double candidateVotes, double topCandidateVotes)
    {
        double votes = support.readVotes();
        return new ContigInfoRow(
                selection, support, candidate,
                share(votes, allVotes), share(votes, candidateVotes), share(votes, topCandidateVotes));
    }

    private enum ContigInfoColumn
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

    private static String formatRanks(Map<ViralContig, Integer> votesRanks, Set<ViralContig> contigs)
    {
        return contigs.stream().map(votesRanks::get).sorted().map(String::valueOf).collect(Collectors.joining(","));
    }

    private static Map<ViralContig, Integer> candidateVotesRanks(List<OncologyGroupRepresentativeSelection> selections)
    {
        return selections.stream()
                .flatMap(selection -> selection.candidates().stream())
                .collect(toMap(RepresentativeContigCandidate::contig, RepresentativeContigCandidate::votesRank));
    }

    // One row per ordered within-oncology-group contig pair: the reads they share, and how many of those fit the subject
    // better by at least each reported margin. Includes support-filtered contigs. Verbose/debug only, for tuning.
    // TODO: write in order: oncology group (lexicographically, ascending), subject rank (ascending), opponent rank (ascending)
    public static void writePairwiseMargins(String file, PairwiseMargins margins, List<OncologyGroupRepresentativeSelection> selections)
    {
        Map<ViralContig, Integer> votesRanks = candidateVotesRanks(selections);

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
                    row.setOrNull(PairwiseMarginsColumn.subject_rank, asString(votesRanks.get(pair.subject())));
                    row.set(PairwiseMarginsColumn.opponent_contig, pair.opponent().name());
                    row.setOrNull(PairwiseMarginsColumn.opponent_rank, asString(votesRanks.get(pair.opponent())));
                    row.set(PairwiseMarginsColumn.shared_reads, margins.sharedReads(pair.subject(), pair.opponent()));
                    for(int margin : REPORTED_MARGINS)
                    {
                        row.set(marginColumn(margin), margins.readsWinningBy(pair.subject(), pair.opponent(), margin));
                    }
                });
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

    private static final List<String> MARGIN_COLUMNS = REPORTED_MARGINS.stream().map(OutputWriter::marginColumn).toList();

    private static String marginColumn(int margin)
    {
        return "reads_m" + margin;
    }

    private static String asString(Object value)
    {
        return value == null ? null : value.toString();
    }
}
