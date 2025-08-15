package com.hartwig.hmftools.panelbuilder;

import static java.lang.String.format;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CANDIDATE_PROBES_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CANDIDATE_REGIONS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_STATS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PANEL_PROBES_FILE_STEM;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.REJECTED_REGIONS_FILE_STEM;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.TARGET_REGIONS_FILE_NAME;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.function.Function;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

// Writes all file output data.
public class OutputWriter implements AutoCloseable
{
    private final DelimFileWriter<Probe> mPanelProbesTsvWriter;
    private final BufferedWriter mPanelProbesBedWriter;
    private final BufferedWriter mPanelProbesFastaWriter;
    private final BufferedWriter mTargetRegionsBedWriter;
    private final DelimFileWriter<RejectedRegion> mRejectedRegionsTsvWriter;
    private final BufferedWriter mRejectedRegionsBedWriter;
    private final BufferedWriter mCandidateRegionsBedWriter;
    @Nullable
    private final DelimFileWriter<Probe> mCandidateProbesTsvWriter;
    @Nullable
    private final ArrayList<Probe> mCandidateProbesBuffer;
    private final DelimFileWriter<Genes.GeneStats> mGeneStatsTsvWriter;

    private static final String TSV_EXT = ".tsv";
    private static final String BED_EXT = ".bed";
    private static final String FASTA_EXT = ".fasta";

    private static final String FLD_SEQUENCE = "Sequence";
    private static final String FLD_QUALITY_SCORE = "QualityScore";
    private static final String FLD_GC_CONTENT = "GCContent";
    private static final String FLD_TARGET_TYPE = "TargetType";
    private static final String FLD_TARGET_EXTRA_INFO = "TargetExtra";
    private static final String FLD_REJECT_REASON = "RejectReason";
    private static final String FLD_EVAL_CRITERIA = "EvalCriteria";
    private static final String FLD_PROBE_COUNT = "ProbeCount";

    private static final List<String> PANEL_PROBES_COLUMNS = List.of(
            FLD_CHROMOSOME, FLD_POSITION_START, FLD_POSITION_END, FLD_SEQUENCE,
            FLD_TARGET_TYPE, FLD_TARGET_EXTRA_INFO,
            FLD_QUALITY_SCORE, FLD_GC_CONTENT);

    private static final List<String> REJECTED_REGIONS_COLUMNS = List.of(
            FLD_CHROMOSOME, FLD_POSITION_START, FLD_POSITION_END,
            FLD_TARGET_TYPE, FLD_TARGET_EXTRA_INFO,
            FLD_REJECT_REASON);

    private static final List<String> CANDIDATE_PROBES_COLUMNS = List.of(
            FLD_CHROMOSOME, FLD_POSITION_START, FLD_POSITION_END, FLD_SEQUENCE,
            FLD_TARGET_TYPE, FLD_TARGET_EXTRA_INFO,
            FLD_QUALITY_SCORE, FLD_GC_CONTENT,
            FLD_EVAL_CRITERIA, FLD_REJECT_REASON);

    private static final List<String> GENE_STATS_COLUMNS = List.of(
            FLD_GENE_NAME,
            FLD_PROBE_COUNT
    );

    private static final int CANDIDATE_PROBES_BUFFER_SIZE = 1_000_000;

    private static final Logger LOGGER = LogManager.getLogger(OutputWriter.class);

    public OutputWriter(final String outputDir, @Nullable final String outputId, boolean verboseOutput) throws IOException
    {
        Function<String, String> outputFilePath = fileName ->
        {
            if(outputId != null)
            {
                fileName = outputId + "." + fileName;
            }
            return Paths.get(outputDir, fileName).toString();
        };

        String panelProbesTsvFile = outputFilePath.apply(PANEL_PROBES_FILE_STEM + TSV_EXT);
        String panelProbesBedFile = outputFilePath.apply(PANEL_PROBES_FILE_STEM + BED_EXT);
        String panelProbesFastaFile = outputFilePath.apply(PANEL_PROBES_FILE_STEM + FASTA_EXT);
        String targetRegionsBedFile = outputFilePath.apply(TARGET_REGIONS_FILE_NAME);
        String rejectedRegionsTsvFile = outputFilePath.apply(REJECTED_REGIONS_FILE_STEM + TSV_EXT);
        String rejectedRegionsBedFile = outputFilePath.apply(REJECTED_REGIONS_FILE_STEM + BED_EXT);
        String candidateRegionsBedFile = outputFilePath.apply(CANDIDATE_REGIONS_FILE_NAME);
        String candidateProbesTsvFile = outputFilePath.apply(CANDIDATE_PROBES_FILE_NAME);
        String geneStatsTsvFile = outputFilePath.apply(GENE_STATS_FILE_NAME);

        mPanelProbesTsvWriter = new DelimFileWriter<>(panelProbesTsvFile, PANEL_PROBES_COLUMNS, OutputWriter::writePanelProbesTsvRow);
        mPanelProbesBedWriter = createBufferedWriter(panelProbesBedFile);
        mPanelProbesFastaWriter = createBufferedWriter(panelProbesFastaFile);

        mTargetRegionsBedWriter = createBufferedWriter(targetRegionsBedFile);

        mRejectedRegionsTsvWriter =
                new DelimFileWriter<>(rejectedRegionsTsvFile, REJECTED_REGIONS_COLUMNS, OutputWriter::writeRejectedRegionsTsvRow);
        mRejectedRegionsBedWriter = createBufferedWriter(rejectedRegionsBedFile);

        mCandidateRegionsBedWriter = createBufferedWriter(candidateRegionsBedFile);

        if(verboseOutput)
        {
            mCandidateProbesTsvWriter =
                    new DelimFileWriter<>(candidateProbesTsvFile, CANDIDATE_PROBES_COLUMNS, OutputWriter::writeCandidateProbesRow);
            mCandidateProbesBuffer = new ArrayList<>(CANDIDATE_PROBES_BUFFER_SIZE);
        }
        else
        {
            mCandidateProbesTsvWriter = null;
            mCandidateProbesBuffer = null;
        }

        mGeneStatsTsvWriter = new DelimFileWriter<>(geneStatsTsvFile, GENE_STATS_COLUMNS, OutputWriter::writeGeneStatsRow);
    }

    public void writePanelProbes(List<Probe> probes) throws IOException
    {
        LOGGER.debug("Writing {} panel probes to file", probes.size());

        // Must be sorted for BED files since some tools expect sorted order.
        probes = probes.stream().sorted(Comparator.comparing(Probe::region, Comparator.nullsLast(Comparator.naturalOrder()))).toList();

        for(Probe probe : probes)
        {
            // A few basic checks that might reveal bugs in the code elsewhere.
            if(!probe.accepted())
            {
                throw new IllegalArgumentException("Should only be writing accepted probes");
            }
            if((probe.region() != null && probe.region().baseLength() != PROBE_LENGTH)
                    || (probe.sequence() != null && probe.sequence().length() != PROBE_LENGTH))
            {
                throw new IllegalArgumentException("Should only be writing probes of length " + PROBE_LENGTH);
            }

            mPanelProbesTsvWriter.writeRow(probe);
            if(probe.region() != null)
            {
                writePanelProbesBedRow(probe);
            }
            writePanelProbesFastaRecord(probe);
        }
    }

    private static void writePanelProbesTsvRow(final Probe probe, DelimFileWriter.Row row)
    {
        ChrBaseRegion region = probe.region();
        row.setOrNull(FLD_CHROMOSOME, region == null ? null : region.chromosome());
        row.setOrNull(FLD_POSITION_START, region == null ? null : region.start());
        row.setOrNull(FLD_POSITION_END, region == null ? null : region.end());
        row.set(FLD_SEQUENCE, probe.sequence());
        row.set(FLD_QUALITY_SCORE, probe.qualityScore());
        row.set(FLD_GC_CONTENT, probe.gcContent());
        row.set(FLD_TARGET_TYPE, probe.metadata().type().name());
        row.set(FLD_TARGET_EXTRA_INFO, probe.metadata().extraInfo());
    }

    private void writePanelProbesBedRow(final Probe probe) throws IOException
    {
        mPanelProbesBedWriter.write(formatBedRow(requireNonNull(probe.region()), probeBedName(probe)));
    }

    private void writePanelProbesFastaRecord(final Probe probe) throws IOException
    {
        String label = getProbeLabel(probe);
        String sequence = probe.sequence();
        if(sequence == null)
        {
            // If this happens there's a code bug.
            throw new IllegalArgumentException("Probe must have sequence data to write FASTA");
        }
        mPanelProbesFastaWriter.write(format(">%s\n%s\n", label, sequence));
    }

    public void writeTargetRegions(List<TargetRegion> regions) throws IOException
    {
        LOGGER.debug("Writing {} target regions to file", regions.size());

        // Must be sorted for BED files since some tools expect sorted order.
        regions = regions.stream().sorted(Comparator.comparing(TargetRegion::region)).toList();

        for(TargetRegion region : regions)
        {
            writeTargetRegionsBedRow(region);
        }
    }

    private void writeTargetRegionsBedRow(final TargetRegion region) throws IOException
    {
        writeTargetRegionBedRow(region, mTargetRegionsBedWriter);
    }

    public void writeRejectedRegions(List<RejectedRegion> regions) throws IOException
    {
        LOGGER.debug("Writing {} rejected regions to file", regions.size());

        // Must be sorted for BED files since some tools expect sorted order.
        regions = regions.stream().sorted(Comparator.comparing(RejectedRegion::region)).toList();

        for(RejectedRegion region : regions)
        {
            mRejectedRegionsTsvWriter.writeRow(region);
            writeRejectedRegionsBedRow(region);
        }
    }

    private static void writeRejectedRegionsTsvRow(final RejectedRegion region, DelimFileWriter.Row row)
    {
        row.set(FLD_CHROMOSOME, region.region().chromosome());
        row.set(FLD_POSITION_START, region.region().start());
        row.set(FLD_POSITION_END, region.region().end());
        row.set(FLD_TARGET_TYPE, region.metadata().type().name());
        row.set(FLD_TARGET_EXTRA_INFO, region.metadata().extraInfo());
        row.set(FLD_REJECT_REASON, region.reason());
    }

    private void writeRejectedRegionsBedRow(final RejectedRegion region) throws IOException
    {
        mRejectedRegionsBedWriter.write(formatBedRow(region.region(), targetMetadataToBedName(region.metadata())));
    }

    public void writeCandidateRegions(List<TargetRegion> regions) throws IOException
    {
        LOGGER.debug("Writing {} candidate regions to file", regions.size());

        // Must be sorted for BED files since some tools expect sorted order.
        regions = regions.stream().sorted(Comparator.comparing(TargetRegion::region)).toList();

        for(TargetRegion region : regions)
        {
            writeTargetRegionBedRow(region, mCandidateRegionsBedWriter);
        }
    }

    public void writeCandidateProbe(final Probe probe)
    {
        if(mCandidateProbesBuffer != null)
        {
            // Buffer probes to improve performance.
            mCandidateProbesBuffer.add(probe);
            checkFlushCandidateProbes(false);
        }
    }

    private void checkFlushCandidateProbes(boolean force)
    {
        if(mCandidateProbesBuffer != null)
        {
            if(mCandidateProbesBuffer.size() >= CANDIDATE_PROBES_BUFFER_SIZE || force)
            {
                writeCandidateProbes(mCandidateProbesBuffer);
                mCandidateProbesBuffer.clear();
            }
        }
    }

    private void writeCandidateProbes(final List<Probe> probes)
    {
        if(mCandidateProbesTsvWriter != null)
        {
            LOGGER.debug("Writing {} candidate probes to file", probes.size());
            probes.forEach(mCandidateProbesTsvWriter::writeRow);
        }
    }

    private static void writeCandidateProbesRow(final Probe probe, DelimFileWriter.Row row)
    {
        ChrBaseRegion region = probe.region();
        row.setOrNull(FLD_CHROMOSOME, region == null ? null : region.chromosome());
        row.setOrNull(FLD_POSITION_START, region == null ? null : region.start());
        row.setOrNull(FLD_POSITION_END, region == null ? null : region.end());
        row.set(FLD_SEQUENCE, probe.sequence());
        row.set(FLD_TARGET_TYPE, probe.metadata().type().name());
        row.set(FLD_TARGET_EXTRA_INFO, probe.metadata().extraInfo());
        row.set(FLD_QUALITY_SCORE, probe.qualityScore());
        row.set(FLD_GC_CONTENT, probe.gcContent());
        row.setOrNull(FLD_EVAL_CRITERIA, requireNonNull(probe.evalCriteria()).toString());
        row.setOrNull(FLD_REJECT_REASON, probe.rejectionReason());
    }

    private static void writeTargetRegionBedRow(final TargetRegion region, BufferedWriter writer) throws IOException
    {
        writer.write(formatBedRow(region.region(), targetMetadataToBedName(region.metadata())));
    }

    private static String probeBedName(final Probe probe)
    {
        double qualityScore = probe.qualityScore();
        double gcContent = probe.gcContent();
        String baseName = targetMetadataToBedName(probe.metadata());
        return format("%s:QS=%.2f:GC=%.2f", baseName, qualityScore, gcContent);
    }

    private static String targetMetadataToBedName(final TargetMetadata info)
    {
        return format("%s:%s", info.type().name(), info.extraInfo());
    }

    private static String formatBedRow(final ChrBaseRegion region, String name)
    {
        return format("%s\t%d\t%d\t%s\n", region.chromosome(), region.start() - 1, region.end(), name);
    }

    private static String getProbeLabel(final Probe probe)
    {
        TargetMetadata metadata = probe.metadata();
        String label = format("%s:%s", metadata.type().name(), metadata.extraInfo());
        if(probe.region() != null)
        {
            label = format("%s:%d", label, probe.region().start());
        }
        return label;
    }

    public void writeGeneStats(final List<Genes.GeneStats> geneStats)
    {
        geneStats.forEach(mGeneStatsTsvWriter::writeRow);
    }

    private static void writeGeneStatsRow(final Genes.GeneStats stats, DelimFileWriter.Row row)
    {
        row.set(FLD_GENE_NAME, stats.geneName());
        row.set(FLD_PROBE_COUNT, stats.probeCount());
    }

    @Override
    public void close() throws IOException
    {
        LOGGER.debug("Flushing and closing output files");

        mPanelProbesTsvWriter.close();
        mPanelProbesBedWriter.close();
        mPanelProbesFastaWriter.close();

        mTargetRegionsBedWriter.close();

        mRejectedRegionsTsvWriter.close();
        mRejectedRegionsBedWriter.close();

        mCandidateRegionsBedWriter.close();

        if(mCandidateProbesTsvWriter != null)
        {
            checkFlushCandidateProbes(true);
            mCandidateProbesTsvWriter.close();
        }

        mGeneStatsTsvWriter.close();
    }
}
