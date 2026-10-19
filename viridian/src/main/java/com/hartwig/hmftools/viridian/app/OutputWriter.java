package com.hartwig.hmftools.viridian.app;

import static java.util.Comparator.comparing;
import static java.util.Comparator.reverseOrder;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toMap;
import static java.util.stream.Collectors.toSet;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.viridian.app.ViridianConfig.CFG_CANDIDATE_READ_SOURCE;
import static com.hartwig.hmftools.viridian.app.ViridianConfig.CFG_EXTRA_VIRUS_DECOY_CONTIGS;
import static com.hartwig.hmftools.viridian.app.ViridianConfig.CFG_KMER_FILTER;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.viridian.detection.DetectedVirus;
import com.hartwig.hmftools.viridian.detection.common.GenomeStats;
import com.hartwig.hmftools.viridian.detection.common.SummaryStats;
import com.hartwig.hmftools.viridian.detection.select.RepresentativeGenome;
import com.hartwig.hmftools.viridian.detection.support.GenomeSupport;
import com.hartwig.hmftools.viridian.integration.IntegrationCandidate;
import com.hartwig.hmftools.viridian.integration.align.ViralInsertAlignment;
import com.hartwig.hmftools.viridian.integration.extract.BreakendSupport;
import com.hartwig.hmftools.viridian.integration.extract.CandidateSv;
import com.hartwig.hmftools.viridian.integration.extract.HostBreakend;
import com.hartwig.hmftools.viridian.integration.extract.InsertRepeat;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.VirusGenome;
import com.hartwig.hmftools.viridian.reporting.VirusReport;

import org.jetbrains.annotations.Nullable;

public class OutputWriter
{
    public static void writeGenomeSupport(
            String file, List<GenomeSupport> genomeSupports, List<RepresentativeGenome> representatives)
    {
        Set<VirusGenome> representativeGenomes =
                representatives.stream().map(RepresentativeGenome::genome).collect(toSet());
        Map<OncologyGroup, GroupVotes> groupVotes = genomeSupports.stream()
                .collect(groupingBy(support -> support.genome().oncologyGroup()))
                .entrySet().stream()
                .collect(toMap(Map.Entry::getKey, entry -> GroupVotes.from(entry.getValue())));
        List<GenomeSupport> rows = genomeSupports.stream()
                .sorted(comparing((GenomeSupport support) -> support.genome().oncologyGroup())
                        .thenComparing(GenomeSupport::readVotes, reverseOrder())
                        .thenComparing(GenomeSupport::genome))
                .toList();

        DelimFileWriter.write(
                file, GENOME_SUPPORT_COLUMNS, rows, (support, row) ->
                {
                    GenomeStats stats = support.stats();
                    VirusGenome genome = support.genome();

                    row.set(GenomeSupportColumn.genome, genome.name());
                    row.set(GenomeSupportColumn.virus_name, genome.virusName());
                    row.set(GenomeSupportColumn.oncology_group, genome.oncologyGroup().name());
                    row.set(GenomeSupportColumn.genome_length, genome.length());
                    row.set(GenomeSupportColumn.read_count, stats.readCount());
                    row.set(GenomeSupportColumn.multi_align_reads, support.multiAlignReads());
                    row.set(GenomeSupportColumn.origin_clipped_reads, stats.originClippedReads());
                    row.set(GenomeSupportColumn.coverage_fraction, stats.coverageFraction());
                    row.set(GenomeSupportColumn.read_votes, support.readVotes());

                    writeSummaryStats(row, DEPTH_STATS_COLUMNS, stats.depth());
                    writeSummaryStats(row, ALIGN_PER_READ_STATS_COLUMNS, support.alignmentsPerRead());
                    writeSummaryStats(row, ALIGNER_SCORE_STATS_COLUMNS, stats.alignerScore());

                    GroupVotes votes = groupVotes.get(genome.oncologyGroup());
                    double readVotes = support.readVotes();
                    row.set(GenomeSupportColumn.filter_status, support.filterStatus().name());
                    row.setOrNull(GenomeSupportColumn.vote_share_pre_filter, share(readVotes, votes.all()));
                    row.setOrNull(GenomeSupportColumn.vote_share_post_filter, share(readVotes, votes.candidates()));
                    row.set(GenomeSupportColumn.is_representative, representativeGenomes.contains(genome));
                });
    }

    private record GroupVotes(
            double all,
            double candidates
    )
    {
        static GroupVotes from(List<GenomeSupport> groupGenomes)
        {
            double all = groupGenomes.stream().mapToDouble(GenomeSupport::readVotes).sum();
            double candidates = groupGenomes.stream()
                    .filter(support -> !support.isFiltered()).mapToDouble(GenomeSupport::readVotes).sum();
            return new GroupVotes(all, candidates);
        }
    }

    private enum GenomeSupportColumn
    {
        genome,
        virus_name,
        oncology_group,
        is_representative,
        filter_status,
        genome_length,
        coverage_fraction,
        read_count,
        read_votes,
        vote_share_pre_filter,
        vote_share_post_filter,
        multi_align_reads,
        origin_clipped_reads
    }

    private static final String DEPTH_STATS_COLUMNS = "depth";
    private static final String ALIGNER_SCORE_STATS_COLUMNS = "aligner_score";

    private static final String ALIGN_PER_READ_STATS_COLUMNS = "align_per_read";

    private static final List<String> GENOME_SUPPORT_COLUMNS =
            Stream.concat(
                            Stream.of(GenomeSupportColumn.values()).map(Enum::name),
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
                    row.set(DetectedVirusColumn.status, detection.representativeStatus().name());

                    GenomeStats stats = detection.representativeStats();
                    VirusGenome genome = stats.genome();
                    row.set(DetectedVirusColumn.genome_name, genome.name());
                    row.set(DetectedVirusColumn.virus_name, genome.virusName());
                    row.set(DetectedVirusColumn.genome_length, genome.length());
                    row.set(DetectedVirusColumn.read_count, stats.readCount());
                    row.set(DetectedVirusColumn.origin_clipped_reads, stats.originClippedReads());
                    row.set(DetectedVirusColumn.coverage_fraction, stats.coverageFraction());

                    writeSummaryStats(row, DEPTH_STATS_COLUMNS, stats.depth());
                    writeSummaryStats(row, ALIGNER_SCORE_STATS_COLUMNS, stats.alignerScore());
                });
    }

    private enum DetectedVirusColumn
    {
        oncology_group,
        virus_name,
        genome_name,
        genome_length,
        coverage_fraction,
        read_count,
        origin_clipped_reads,
        status,
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

    public static void writeIntegrationCandidates(String file, List<IntegrationCandidate> candidates)
    {
        DelimFileWriter.write(
                file, IntegrationCandidateColumn.values(), candidates, (candidate, row) ->
                {
                    row.set(IntegrationCandidateColumn.is_integration, candidate.isIntegration());

                    CandidateSv sv = candidate.sv();
                    row.set(IntegrationCandidateColumn.sv_type, sv.type().name());
                    row.set(IntegrationCandidateColumn.filter, sv.filter());

                    HostBreakend startBreakend = sv.startBreakend();
                    BreakendSupport startSupport = startBreakend.support();
                    row.set(IntegrationCandidateColumn.start_id, startBreakend.id());
                    row.set(IntegrationCandidateColumn.start_chromosome, startBreakend.basePosition().Chromosome);
                    row.set(IntegrationCandidateColumn.start_position, startBreakend.basePosition().Position);
                    row.set(IntegrationCandidateColumn.start_orientation, startBreakend.orientation().asByte());
                    row.set(IntegrationCandidateColumn.start_tumor_variant_frags, startSupport.tumorVariantFragments());
                    row.set(IntegrationCandidateColumn.start_tumor_ref_frags, startSupport.tumorReferenceFragments());
                    row.set(IntegrationCandidateColumn.start_tumor_af, startSupport.tumorAlleleFrequency());
                    row.setOrNull(
                            IntegrationCandidateColumn.start_normal_variant_frags.name(),
                            startSupport.normalVariantFragments());
                    row.setOrNull(
                            IntegrationCandidateColumn.start_normal_ref_frags.name(),
                            startSupport.normalReferenceFragments());

                    // Note unset columns are written as null, leaving a SGL's second breakend blank.
                    HostBreakend endBreakend = sv.endBreakend();
                    if(endBreakend != null)
                    {
                        BreakendSupport endSupport = endBreakend.support();
                        row.set(IntegrationCandidateColumn.end_id, endBreakend.id());
                        row.set(IntegrationCandidateColumn.end_chromosome, endBreakend.basePosition().Chromosome);
                        row.set(IntegrationCandidateColumn.end_position, endBreakend.basePosition().Position);
                        row.set(IntegrationCandidateColumn.end_orientation, endBreakend.orientation().asByte());
                        row.set(IntegrationCandidateColumn.end_tumor_variant_frags, endSupport.tumorVariantFragments());
                        row.set(IntegrationCandidateColumn.end_tumor_ref_frags, endSupport.tumorReferenceFragments());
                        row.set(IntegrationCandidateColumn.end_tumor_af, endSupport.tumorAlleleFrequency());
                        row.setOrNull(
                                IntegrationCandidateColumn.end_normal_variant_frags.name(),
                                endSupport.normalVariantFragments());
                        row.setOrNull(
                                IntegrationCandidateColumn.end_normal_ref_frags.name(),
                                endSupport.normalReferenceFragments());
                    }

                    row.set(IntegrationCandidateColumn.line_site, sv.lineSite());
                    InsertRepeat insertRepeat = sv.insertRepeat();
                    if(insertRepeat != null)
                    {
                        row.set(IntegrationCandidateColumn.repeat_class, insertRepeat.repeatClass());
                        row.set(IntegrationCandidateColumn.repeat_type, insertRepeat.repeatType());
                        row.set(IntegrationCandidateColumn.repeat_coverage, insertRepeat.coverage());
                    }
                    row.set(IntegrationCandidateColumn.insert_host_alignments, sv.insertHostAlignments());
                    row.set(IntegrationCandidateColumn.insert_seq_length, sv.insertSequence().length());
                    row.set(IntegrationCandidateColumn.insert_sequence, sv.insertSequence());

                    ViralInsertAlignment alignment = candidate.alignment();
                    if(alignment != null)
                    {
                        VirusGenome genome = alignment.genome();
                        row.set(IntegrationCandidateColumn.virus_name, genome.virusName());
                        row.set(IntegrationCandidateColumn.oncology_group, genome.oncologyGroup().name());
                        row.set(IntegrationCandidateColumn.virus_contig, genome.name());
                        row.set(IntegrationCandidateColumn.virus_position, alignment.position());
                        row.set(IntegrationCandidateColumn.virus_orientation, alignment.orientation().asByte());
                        row.set(IntegrationCandidateColumn.align_cigar, alignment.cigar().toString());
                        row.set(IntegrationCandidateColumn.aligned_length, alignment.alignedLength());
                        row.set(IntegrationCandidateColumn.aligner_score, alignment.alignerScore());
                        row.set(IntegrationCandidateColumn.score_per_aligned_base, alignment.scorePerAlignedBase());
                        row.set(IntegrationCandidateColumn.aligned_edit_distance, alignment.alignedEditDistance());
                    }
                });
    }

    private enum IntegrationCandidateColumn
    {
        is_integration,
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
        aligned_edit_distance
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
                    row.set(VirusReportColumn.status, report.status().name());
                    row.set(VirusReportColumn.integrations, report.integrations());
                    row.set(VirusReportColumn.detected, report.isDetected());

                    DetectedVirus detection = report.detection();
                    if(detection != null)
                    {
                        GenomeStats stats = detection.representativeStats();
                        VirusGenome genome = stats.genome();
                        row.set(VirusReportColumn.representative_status, detection.representativeStatus().name());
                        row.set(VirusReportColumn.representative_genome, genome.name());
                        row.set(VirusReportColumn.virus_name, genome.virusName());
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
        status,
        integrations,
        detected,
        representative_status,
        representative_genome,
        virus_name,
        coverage_fraction,
        mean_depth,
        copies_per_tumor_cell
    }

    private static final List<String> VIRUS_REPORT_COLUMNS =
            Stream.of(VirusReportColumn.values()).map(Enum::name).toList();

    public static void writeRunManifest(String file, ViridianConfig config)
    {
        List<String> columns = List.of(CFG_CANDIDATE_READ_SOURCE, CFG_EXTRA_VIRUS_DECOY_CONTIGS, CFG_KMER_FILTER);
        DelimFileWriter.write(
                file, columns, List.of(config), (runConfig, row) ->
                {
                    row.set(CFG_CANDIDATE_READ_SOURCE, runConfig.candidateReadSource().name());
                    row.set(
                            CFG_EXTRA_VIRUS_DECOY_CONTIGS,
                            String.join(ITEM_DELIM, runConfig.extraVirusDecoyContigs().stream().sorted().toList()));
                    row.set(CFG_KMER_FILTER, runConfig.kmerFilterEnabled());
                });
    }
}
