package com.hartwig.hmftools.virusdetect;

import static java.util.Comparator.comparing;
import static java.util.Comparator.comparingInt;
import static java.util.stream.Collectors.toMap;

import static com.hartwig.hmftools.virusdetect.VirusConstants.CHALLENGE_MARGIN_SWEEP;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class VirusOutputWriter
{
    private static final Logger LOGGER = LogManager.getLogger(VirusOutputWriter.class);

    public static void writeContigStats(
            String file, Collection<ContigStats> stats, List<ContigClassification> classifications)
    {
        Map<ViralContig, ContigClassification> byContig = classifications.stream()
                .collect(toMap(ContigClassification::contig, classification -> classification));

        List<ContigStats> ordered = stats.stream()
                .sorted(comparingInt(ContigStats::readCount).reversed().thenComparing(stat -> stat.contig().name()))
                .toList();

        DelimFileWriter.write(
                file, Column.values(), ordered, (stat, row) ->
                {
                    ViralContig contig = stat.contig();
                    row.set(Column.contig, contig.name());
                    row.set(Column.virus_name, contig.virusName());
                    row.set(Column.oncology_group, contig.oncologyGroup());
                    row.set(Column.contig_length, contig.length());
                    row.set(Column.read_count, stat.readCount());
                    row.set(Column.multi_align_reads, stat.multiAlignReads());
                    row.set(Column.origin_clipped_reads, stat.originClippedReads());
                    row.set(Column.coverage_fraction, stat.coverageFraction());
                    row.set(Column.read_votes, stat.readVotes());

                    SummaryStats alignPerRead = stat.alignPerRead();
                    row.set(Column.align_per_read_mean, alignPerRead.mean());
                    row.set(Column.align_per_read_min, alignPerRead.min());
                    row.set(Column.align_per_read_p5, alignPerRead.p5());
                    row.set(Column.align_per_read_p25, alignPerRead.p25());
                    row.set(Column.align_per_read_p50, alignPerRead.p50());
                    row.set(Column.align_per_read_p75, alignPerRead.p75());
                    row.set(Column.align_per_read_p95, alignPerRead.p95());
                    row.set(Column.align_per_read_max, alignPerRead.max());

                    SummaryStats depth = stat.depth();
                    row.set(Column.depth_mean, depth.mean());
                    row.set(Column.depth_min, depth.min());
                    row.set(Column.depth_p5, depth.p5());
                    row.set(Column.depth_p25, depth.p25());
                    row.set(Column.depth_p50, depth.p50());
                    row.set(Column.depth_p75, depth.p75());
                    row.set(Column.depth_p95, depth.p95());
                    row.set(Column.depth_max, depth.max());

                    SummaryStats alignerScore = stat.alignerScore();
                    row.set(Column.aligner_score_mean, alignerScore.mean());
                    row.set(Column.aligner_score_min, alignerScore.min());
                    row.set(Column.aligner_score_p5, alignerScore.p5());
                    row.set(Column.aligner_score_p25, alignerScore.p25());
                    row.set(Column.aligner_score_p50, alignerScore.p50());
                    row.set(Column.aligner_score_p75, alignerScore.p75());
                    row.set(Column.aligner_score_p95, alignerScore.p95());
                    row.set(Column.aligner_score_max, alignerScore.max());

                    ContigClassification classification = byContig.get(stat.contig());
                    row.set(Column.filter_status, classification.filterStatus().name());
                    row.setOrNull(Column.votes_rank, asString(classification.votesRank()));
                    row.setOrNull(Column.vote_share_pre_filter, classification.preFilterVoteShare());
                    row.setOrNull(Column.vote_share_post_filter, classification.postFilterVoteShare());
                    row.setOrNull(Column.vote_share_ratio, classification.voteShareRatio());
                    row.setOrNull(Column.comparable, asString(classification.comparable()));
                    row.setOrNull(Column.challenges_ranks, candidateRanks(classification, classification.challengesRanks()));
                    row.setOrNull(Column.challenged_by_ranks, candidateRanks(classification, classification.challengedByRanks()));
                    row.setOrNull(Column.role, asString(classification.role()));
                    row.setOrNull(Column.oncology_group_outcome, asString(classification.oncologyGroupOutcome()));
                    row.setOrNull(Column.oncology_group_sub_outcome, asString(classification.oncologyGroupSubOutcome()));
                });

        LOGGER.info("wrote {} contig stats to {}", ordered.size(), file);
    }

    // One row per ordered within-oncology-group contig pair, how decisively the subject fits shared reads better than the
    // opponent as a challenge share at each considered margin value. Includes filtered contigs. Verbose/debug only, for tuning.
    public static void writePairwiseMargins(
            String file, PairwiseMargins pairwise, RepresentativeSelectionResult selection)
    {
        Map<ViralContig, ContigClassification> byContig = selection.classifications().stream()
                .collect(toMap(ContigClassification::contig, classification -> classification));

        List<PairwiseMargins.ContigPair> pairs = pairwise.pairs().stream()
                .sorted(comparing((PairwiseMargins.ContigPair pair) -> pair.subject().name())
                        .thenComparing(pair -> pair.opponent().name()))
                .toList();

        List<String> columns = Stream.concat(
                        Stream.of(PairwiseColumn.values()).map(Enum::name),
                        marginColumns().stream())
                .toList();

        DelimFileWriter.write(
                file, columns, pairs, (pair, row) ->
                {
                    String oncologyGroup = pair.subject().oncologyGroup();
                    double voteTotal = selection.oncologyGroupVoteTotals().getOrDefault(oncologyGroup, 0.0);
                    row.set(PairwiseColumn.oncology_group, oncologyGroup);
                    row.set(PairwiseColumn.subject_contig, pair.subject().name());
                    row.setOrNull(PairwiseColumn.subject_rank, asString(votesRank(byContig, pair.subject())));
                    row.set(PairwiseColumn.opponent_contig, pair.opponent().name());
                    row.setOrNull(PairwiseColumn.opponent_rank, asString(votesRank(byContig, pair.opponent())));
                    row.set(PairwiseColumn.shared_reads, pairwise.sharedReads(pair.subject(), pair.opponent()));
                    for(int margin : CHALLENGE_MARGIN_SWEEP)
                    {
                        double challengeShare =
                                voteTotal > 0 ? pairwise.challengeReads(pair.subject(), pair.opponent(), margin) / voteTotal : 0.0;
                        row.set(marginColumn(margin), challengeShare);
                    }
                });

        LOGGER.info("wrote {} pairwise margin rows to {}", pairs.size(), file);
    }

    private static Integer votesRank(Map<ViralContig, ContigClassification> byContig, ViralContig contig)
    {
        ContigClassification classification = byContig.get(contig);
        return classification == null ? null : classification.votesRank();
    }

    private static List<String> marginColumns()
    {
        List<String> columns = new ArrayList<>();
        for(int margin : CHALLENGE_MARGIN_SWEEP)
        {
            columns.add(marginColumn(margin));
        }
        return columns;
    }

    private static String marginColumn(int margin)
    {
        return "share_m" + margin;
    }

    // Rank lists apply only to candidates (which carry a role); a non-candidate has no challenge relations to report.
    private static String candidateRanks(ContigClassification classification, List<Integer> ranks)
    {
        if(classification.role() == null)
        {
            return null;
        }
        return ranks.stream().map(String::valueOf).collect(Collectors.joining(","));
    }

    private static String asString(Object value)
    {
        return value == null ? null : value.toString();
    }

    private enum Column
    {
        contig,
        virus_name,
        oncology_group,
        oncology_group_outcome,
        oncology_group_sub_outcome,
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

    private enum PairwiseColumn
    {
        oncology_group,
        subject_contig,
        subject_rank,
        opponent_contig,
        opponent_rank,
        shared_reads
    }
}
