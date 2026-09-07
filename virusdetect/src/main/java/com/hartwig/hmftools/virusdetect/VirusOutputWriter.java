package com.hartwig.hmftools.virusdetect;

import static java.util.Comparator.comparingInt;

import java.util.Collection;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class VirusOutputWriter
{
    private static final Logger LOGGER = LogManager.getLogger(VirusOutputWriter.class);

    public static void writeContigStats(String file, Collection<ContigStats> stats, ViralReference reference)
    {
        List<ContigStats> ordered = stats.stream()
                .sorted(comparingInt(ContigStats::readCount).reversed().thenComparing(ContigStats::contig))
                .toList();

        DelimFileWriter.write(
                file, Column.values(), ordered, (stat, row) ->
                {
                    ViralContig contig = reference.contig(stat.contig());
                    row.set(Column.contig, stat.contig());
                    row.set(Column.virus_name, contig.virusName());
                    row.set(Column.oncology_group, contig.oncologyGroup());
                    row.set(Column.contig_length, stat.contigLength());
                    row.set(Column.read_count, stat.readCount());
                    row.set(Column.covered_bases, stat.coveredBases());
                    row.set(Column.coverage_fraction, stat.coverageFraction());

                    SummaryStats depth = stat.depth();
                    row.set(Column.depth_mean, depth.mean());
                    row.set(Column.depth_min, depth.min());
                    row.set(Column.depth_p5, depth.p5());
                    row.set(Column.depth_p50, depth.p50());
                    row.set(Column.depth_p95, depth.p95());
                    row.set(Column.depth_max, depth.max());

                    row.set(Column.mean_aligner_score, stat.meanAlignerScore());
                    row.set(Column.read_votes, stat.readVotes());
                    row.set(Column.reads_best_in_rivals, stat.readsBestInRivals());

                    Optional<SummaryStats> margins = stat.margins();
                    row.setOrNull(Column.margin_mean, margins.map(SummaryStats::mean).orElse(null));
                    row.setOrNull(Column.margin_min, margins.map(SummaryStats::min).orElse(null));
                    row.setOrNull(Column.margin_p5, margins.map(SummaryStats::p5).orElse(null));
                    row.setOrNull(Column.margin_p50, margins.map(SummaryStats::p50).orElse(null));
                    row.setOrNull(Column.margin_p95, margins.map(SummaryStats::p95).orElse(null));
                    row.setOrNull(Column.margin_max, margins.map(SummaryStats::max).orElse(null));
                });

        LOGGER.info("wrote {} contig stats to {}", ordered.size(), file);
    }

    private enum Column
    {
        contig,
        virus_name,
        oncology_group,
        contig_length,
        read_count,
        covered_bases,
        coverage_fraction,
        depth_mean,
        depth_min,
        depth_p5,
        depth_p50,
        depth_p95,
        depth_max,
        mean_aligner_score,
        read_votes,
        reads_best_in_rivals,
        margin_mean,
        margin_min,
        margin_p5,
        margin_p50,
        margin_p95,
        margin_max
    }
}
