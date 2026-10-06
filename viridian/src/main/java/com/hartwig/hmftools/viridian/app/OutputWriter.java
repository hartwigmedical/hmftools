package com.hartwig.hmftools.viridian.app;

import static java.util.Comparator.comparing;
import static java.util.Comparator.naturalOrder;
import static java.util.Comparator.nullsLast;
import static java.util.Comparator.reverseOrder;
import static java.util.stream.Collectors.toMap;

import static com.hartwig.hmftools.viridian.common.ViridianConstants.REPORTED_MARGINS;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.viridian.detection.DetectedVirus;
import com.hartwig.hmftools.viridian.detection.common.ContigStats;
import com.hartwig.hmftools.viridian.detection.common.SummaryStats;
import com.hartwig.hmftools.viridian.detection.select.OncologyGroupRepresentativeSelection;
import com.hartwig.hmftools.viridian.detection.select.PairwiseMargins;
import com.hartwig.hmftools.viridian.detection.select.RepresentativeContigCandidate;
import com.hartwig.hmftools.viridian.detection.support.ContigSupport;
import com.hartwig.hmftools.viridian.integration.Integration;
import com.hartwig.hmftools.viridian.integration.align.ViralInsertAlignment;
import com.hartwig.hmftools.viridian.integration.extract.BreakendSupport;
import com.hartwig.hmftools.viridian.integration.extract.HostBreakend;
import com.hartwig.hmftools.viridian.integration.extract.HostVariantCandidate;
import com.hartwig.hmftools.viridian.integration.extract.InsertRepeat;
import com.hartwig.hmftools.viridian.reference.ViralContig;
import com.hartwig.hmftools.viridian.reporting.VirusReport;

import org.jetbrains.annotations.Nullable;

public class OutputWriter
{
    public static void writeContigSupport(String file, List<OncologyGroupRepresentativeSelection> selections)
    {
        Map<ViralContig, Integer> votesRanks = candidateVotesRanks(selections);
        List<ContigSupportRow> rows = selections.stream()
                .flatMap(OutputWriter::contigSupportRows)
                .sorted(comparing((ContigSupportRow row) -> row.support().contig().oncologyGroup())
                        .thenComparing(row -> row.support().readVotes(), reverseOrder())
                        .thenComparing(row -> row.support().contig()))
                .toList();

        DelimFileWriter.write(
                file, CONTIG_SUPPORT_COLUMNS, rows, (contigRow, row) ->
                {
                    ContigSupport support = contigRow.support();
                    ContigStats stats = support.stats();
                    ViralContig contig = support.contig();

                    row.set(ContigSupportColumn.contig, contig.name());
                    row.set(ContigSupportColumn.virus_name, contig.virusName());
                    row.set(ContigSupportColumn.oncology_group, contig.oncologyGroup().name());
                    row.set(ContigSupportColumn.contig_length, contig.length());
                    row.set(ContigSupportColumn.read_count, stats.readCount());
                    row.set(ContigSupportColumn.multi_align_reads, support.multiAlignReads());
                    row.set(ContigSupportColumn.origin_clipped_reads, stats.originClippedReads());
                    row.set(ContigSupportColumn.coverage_fraction, stats.coverageFraction());
                    row.set(ContigSupportColumn.read_votes, support.readVotes());

                    writeSummaryStats(row, DEPTH_STATS_COLUMNS, stats.depth());
                    writeSummaryStats(row, ALIGN_PER_READ_STATS_COLUMNS, support.alignmentsPerRead());
                    writeSummaryStats(row, ALIGNER_SCORE_STATS_COLUMNS, stats.alignerScore());

                    row.set(ContigSupportColumn.oncology_group_resolution, contigRow.selection().resolution().name());
                    row.set(ContigSupportColumn.oncology_group_outcome, contigRow.selection().outcome().name());
                    row.set(ContigSupportColumn.filter_status, support.filterStatus().name());
                    row.setOrNull(ContigSupportColumn.vote_share_pre_filter, contigRow.preFilterVoteShare());
                    row.setOrNull(ContigSupportColumn.vote_share_post_filter, contigRow.postFilterVoteShare());
                    row.setOrNull(ContigSupportColumn.vote_share_ratio, contigRow.voteShareOfTop());

                    // Note unset columns are written as null.
                    RepresentativeContigCandidate candidate = contigRow.candidate();
                    if(candidate != null)
                    {
                        row.set(ContigSupportColumn.votes_rank, candidate.votesRank());
                        row.set(ContigSupportColumn.comparable, candidate.comparable());
                        row.set(ContigSupportColumn.role, candidate.role().name());
                        row.set(ContigSupportColumn.challenges_ranks, formatRanks(votesRanks, candidate.challenges()));
                        row.set(ContigSupportColumn.challenged_by_ranks, formatRanks(votesRanks, candidate.challengedBy()));
                    }
                });
    }

    private record ContigSupportRow(
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

    private static Stream<ContigSupportRow> contigSupportRows(OncologyGroupRepresentativeSelection selection)
    {
        double candidateVotes = selection.candidates().stream().mapToDouble(c -> c.support().readVotes()).sum();
        double rejectedVotes = selection.rejected().stream().mapToDouble(ContigSupport::readVotes).sum();
        double topCandidateVotes = selection.candidates().stream()
                .mapToDouble(candidate -> candidate.support().readVotes()).max().orElse(0.0);

        Stream<ContigSupportRow> candidates = selection.candidates().stream()
                .map(candidate -> contigSupportRow(
                        selection, candidate.support(), candidate, candidateVotes + rejectedVotes,
                        candidateVotes, topCandidateVotes));
        Stream<ContigSupportRow> rejected = selection.rejected().stream()
                .map(support -> contigSupportRow(
                        selection, support, null, candidateVotes + rejectedVotes,
                        candidateVotes, topCandidateVotes));

        return Stream.concat(candidates, rejected);
    }

    private static ContigSupportRow contigSupportRow(
            OncologyGroupRepresentativeSelection selection, ContigSupport support, @Nullable RepresentativeContigCandidate candidate,
            double allVotes, double candidateVotes, double topCandidateVotes)
    {
        double votes = support.readVotes();
        return new ContigSupportRow(
                selection, support, candidate,
                share(votes, allVotes), share(votes, candidateVotes), share(votes, topCandidateVotes));
    }

    private enum ContigSupportColumn
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

    private static final String DEPTH_STATS_COLUMNS = "depth";
    private static final String ALIGNER_SCORE_STATS_COLUMNS = "aligner_score";

    private static final String ALIGN_PER_READ_STATS_COLUMNS = "align_per_read";

    private static final List<String> CONTIG_SUPPORT_COLUMNS =
            Stream.concat(
                            Stream.of(ContigSupportColumn.values()).map(Enum::name),
                            Stream.of(DEPTH_STATS_COLUMNS, ALIGN_PER_READ_STATS_COLUMNS, ALIGNER_SCORE_STATS_COLUMNS)
                                    .flatMap(group -> SummaryStats.FIELD_NAMES.stream()
                                            .map(field -> summaryStatsColumn(group, field))))
                    .toList();

    public static void writeDetectedViruses(String file, List<DetectedVirus> detections)
    {
        List<DetectedVirus> rows = detections.stream().sorted(comparing(DetectedVirus::oncologyGroup)).toList();

        DelimFileWriter.write(
                file, VIRUS_DETECTION_COLUMNS, rows, (detection, row) ->
                {
                    row.set(DetectedVirusColumn.oncology_group, detection.oncologyGroup().name());
                    row.set(DetectedVirusColumn.oncology_group_resolution, detection.resolution().name());
                    row.set(DetectedVirusColumn.oncology_group_outcome, detection.outcome().name());
                    row.set(DetectedVirusColumn.group_read_count, detection.groupReadCount());
                    row.set(DetectedVirusColumn.aligned_contig_count, detection.alignedContigCount());
                    row.set(DetectedVirusColumn.candidate_count, detection.candidateCount());
                    row.set(DetectedVirusColumn.comparable_candidate_count, detection.comparableCandidateCount());

                    // Note unset columns are written as null, leaving an unmeasured group's columns blank.
                    ContigStats stats = detection.representativeContigStats();
                    if(stats != null)
                    {
                        ViralContig contig = stats.contig();
                        row.set(DetectedVirusColumn.representative_contig, contig.name());
                        row.set(DetectedVirusColumn.virus_name, contig.virusName());
                        row.set(DetectedVirusColumn.contig_length, contig.length());
                        row.set(DetectedVirusColumn.read_count, stats.readCount());
                        row.set(DetectedVirusColumn.origin_clipped_reads, stats.originClippedReads());
                        row.set(DetectedVirusColumn.coverage_fraction, stats.coverageFraction());

                        writeSummaryStats(row, DEPTH_STATS_COLUMNS, stats.depth());
                        writeSummaryStats(row, ALIGNER_SCORE_STATS_COLUMNS, stats.alignerScore());
                    }
                });
    }

    private enum DetectedVirusColumn
    {
        oncology_group,
        oncology_group_resolution,
        oncology_group_outcome,
        group_read_count,
        aligned_contig_count,
        candidate_count,
        comparable_candidate_count,
        representative_contig,
        virus_name,
        contig_length,
        coverage_fraction,
        read_count,
        origin_clipped_reads
    }

    private static final List<String> VIRUS_DETECTION_COLUMNS =
            Stream.concat(
                            Stream.of(DetectedVirusColumn.values()).map(Enum::name),
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
        return group + "_" + field;
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
    public static void writePairwiseMargins(String file, PairwiseMargins margins, List<OncologyGroupRepresentativeSelection> selections)
    {
        Map<ViralContig, Integer> votesRanks = candidateVotesRanks(selections);

        // A pair's contigs share an oncology group, and only its candidates are ranked, so unranked contigs sort last.
        List<PairwiseMargins.ContigPair> pairs = margins.pairs().stream()
                .sorted(comparing((PairwiseMargins.ContigPair pair) -> pair.subject().oncologyGroup())
                        .thenComparing(pair -> votesRanks.get(pair.subject()), nullsLast(naturalOrder()))
                        .thenComparing(PairwiseMargins.ContigPair::subject)
                        .thenComparing(pair -> votesRanks.get(pair.opponent()), nullsLast(naturalOrder()))
                        .thenComparing(PairwiseMargins.ContigPair::opponent))
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

    public static void writeIntegrations(String file, List<Integration> integrations)
    {
        DelimFileWriter.write(
                file, IntegrationColumn.values(), integrations, (integration, row) ->
                {
                    HostVariantCandidate candidate = integration.hostVariant();
                    row.set(IntegrationColumn.sv_type, candidate.type().name());
                    row.set(IntegrationColumn.filter, candidate.filter());

                    HostBreakend startBreakend = candidate.startBreakend();
                    BreakendSupport startSupport = startBreakend.support();
                    row.set(IntegrationColumn.start_id, startBreakend.id());
                    row.set(IntegrationColumn.start_chromosome, startBreakend.basePosition().Chromosome);
                    row.set(IntegrationColumn.start_position, startBreakend.basePosition().Position);
                    row.set(IntegrationColumn.start_orientation, startBreakend.orientation().asByte());
                    row.set(IntegrationColumn.start_tumor_variant_frags, startSupport.tumorVariantFragments());
                    row.set(IntegrationColumn.start_tumor_ref_frags, startSupport.tumorReferenceFragments());
                    row.set(IntegrationColumn.start_tumor_af, startSupport.tumorAlleleFrequency());
                    row.setOrNull(IntegrationColumn.start_normal_variant_frags.name(), startSupport.normalVariantFragments());
                    row.setOrNull(IntegrationColumn.start_normal_ref_frags.name(), startSupport.normalReferenceFragments());

                    // Note unset columns are written as null, leaving a SGL's second breakend blank.
                    HostBreakend endBreakend = candidate.endBreakend();
                    if(endBreakend != null)
                    {
                        BreakendSupport endSupport = endBreakend.support();
                        row.set(IntegrationColumn.end_id, endBreakend.id());
                        row.set(IntegrationColumn.end_chromosome, endBreakend.basePosition().Chromosome);
                        row.set(IntegrationColumn.end_position, endBreakend.basePosition().Position);
                        row.set(IntegrationColumn.end_orientation, endBreakend.orientation().asByte());
                        row.set(IntegrationColumn.end_tumor_variant_frags, endSupport.tumorVariantFragments());
                        row.set(IntegrationColumn.end_tumor_ref_frags, endSupport.tumorReferenceFragments());
                        row.set(IntegrationColumn.end_tumor_af, endSupport.tumorAlleleFrequency());
                        row.setOrNull(IntegrationColumn.end_normal_variant_frags.name(), endSupport.normalVariantFragments());
                        row.setOrNull(IntegrationColumn.end_normal_ref_frags.name(), endSupport.normalReferenceFragments());
                    }

                    row.set(IntegrationColumn.line_site, candidate.lineSite());
                    InsertRepeat insertRepeat = candidate.insertRepeat();
                    if(insertRepeat != null)
                    {
                        row.set(IntegrationColumn.repeat_class, insertRepeat.repeatClass());
                        row.set(IntegrationColumn.repeat_type, insertRepeat.repeatType());
                        row.set(IntegrationColumn.repeat_coverage, insertRepeat.coverage());
                    }
                    row.set(IntegrationColumn.insert_host_alignments, candidate.insertHostAlignments());
                    row.set(IntegrationColumn.insert_seq_length, candidate.insertSequence().length());
                    row.set(IntegrationColumn.insert_sequence, candidate.insertSequence());

                    ViralInsertAlignment alignment = integration.alignment();
                    if(alignment != null)
                    {
                        ViralContig contig = alignment.contig();
                        row.set(IntegrationColumn.virus_name, contig.virusName());
                        row.set(IntegrationColumn.oncology_group, contig.oncologyGroup().name());
                        row.set(IntegrationColumn.virus_contig, contig.name());
                        row.set(IntegrationColumn.virus_position, alignment.position());
                        row.set(IntegrationColumn.virus_orientation, alignment.orientation().asByte());
                        row.set(IntegrationColumn.align_cigar, alignment.cigar().toString());
                        row.set(IntegrationColumn.aligned_length, alignment.alignedLength());
                        row.set(IntegrationColumn.aligner_score, alignment.alignerScore());
                        row.set(IntegrationColumn.score_per_aligned_base, alignment.scorePerAlignedBase());
                        row.set(IntegrationColumn.aligned_edit_distance, alignment.alignedEditDistance());
                        row.set(IntegrationColumn.plausible, integration.isPlausible());
                    }
                });
    }

    private enum IntegrationColumn
    {
        sv_type,
        filter,
        start_id,
        start_chromosome,
        start_position,
        start_orientation,
        end_id,
        end_chromosome,
        end_position,
        end_orientation,
        start_tumor_variant_frags,
        start_tumor_ref_frags,
        start_tumor_af,
        start_normal_variant_frags,
        start_normal_ref_frags,
        end_tumor_variant_frags,
        end_tumor_ref_frags,
        end_tumor_af,
        end_normal_variant_frags,
        end_normal_ref_frags,
        line_site,
        repeat_class,
        repeat_type,
        repeat_coverage,
        insert_host_alignments,
        insert_seq_length,
        insert_sequence,
        virus_name,
        oncology_group,
        virus_contig,
        virus_position,
        virus_orientation,
        align_cigar,
        aligned_length,
        aligner_score,
        score_per_aligned_base,
        aligned_edit_distance,
        plausible
    }

    public static void writeVirusReports(String file, List<VirusReport> reports)
    {
        List<VirusReport> rows = reports.stream()
                .filter(report -> report.reportingType() != null)
                .sorted(comparing(VirusReport::oncologyGroup))
                .toList();

        DelimFileWriter.write(
                file, VIRUS_REPORT_COLUMNS, rows, (report, row) ->
                {
                    row.set(VirusReportColumn.oncology_group, report.oncologyGroup().name());
                    row.set(VirusReportColumn.reporting_type, report.reportingType().toString());
                    row.set(VirusReportColumn.driver_likelihood, report.driverLikelihood().toString());
                    row.set(VirusReportColumn.reported, report.isReported());
                    row.set(VirusReportColumn.reason, report.reason().name());
                    row.set(VirusReportColumn.integrations, report.integrations());
                    row.set(VirusReportColumn.present, report.isPresent());

                    ContigStats stats = report.representativeStats();
                    if(stats != null)
                    {
                        ViralContig contig = stats.contig();
                        row.set(VirusReportColumn.representative_contig, contig.name());
                        row.set(VirusReportColumn.virus_name, contig.virusName());
                        row.set(VirusReportColumn.coverage_fraction, stats.coverageFraction());
                        row.set(VirusReportColumn.mean_depth, stats.depth().mean());
                    }
                    row.setOrNull(VirusReportColumn.copies_per_tumor_cell, report.copiesPerTumorCell());
                });
    }

    private enum VirusReportColumn
    {
        oncology_group,
        reporting_type,
        driver_likelihood,
        reported,
        reason,
        integrations,
        present,
        representative_contig,
        virus_name,
        coverage_fraction,
        mean_depth,
        copies_per_tumor_cell
    }

    private static final List<String> VIRUS_REPORT_COLUMNS =
            Stream.of(VirusReportColumn.values()).map(Enum::name).toList();

    private static String asString(@Nullable Object value)
    {
        return value == null ? null : value.toString();
    }
}
