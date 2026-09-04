package com.hartwig.hmftools.virusdetect;

import static java.util.Comparator.comparingInt;

import java.util.Collection;
import java.util.List;

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
                    row.set(Column.min_depth, stat.minDepth());
                    row.set(Column.max_depth, stat.maxDepth());
                    row.set(Column.mean_depth, stat.meanDepth());
                    row.set(Column.mean_aligner_score, stat.meanAlignerScore());
                    row.set(Column.read_votes, stat.readVotes());
                    row.set(Column.reads_best_in_rivals, stat.readsBestInRivals());

                    // A contig no read contests holds no margins; write null rather than a misleading zero.
                    boolean contested = stat.readsBestInRivals() > 0;
                    row.setOrNull(Column.margin_mean, contested ? stat.marginMean() : null);
                    row.setOrNull(Column.margin_median, contested ? stat.marginMedian() : null);
                    row.setOrNull(Column.margin_p90, contested ? stat.marginP90() : null);
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
        min_depth,
        max_depth,
        mean_depth,
        mean_aligner_score,
        read_votes,
        reads_best_in_rivals,
        margin_mean,
        margin_median,
        margin_p90
    }
}
