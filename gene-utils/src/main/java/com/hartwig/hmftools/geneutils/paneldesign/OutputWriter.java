package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

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
    private final BufferedWriter mTargetRegionsWriter;
    private final DelimFileWriter<RejectedRegion> mRejectedRegionsTsvWriter;
    private final BufferedWriter mRejectedRegionsBedWriter;
    @Nullable
    private final DelimFileWriter<Probe> mCandidateProbesWriter;
    @Nullable
    private final List<Probe> mCandidateProbesBuffer;

    private static final String TSV_EXT = ".tsv";
    private static final String BED_EXT = ".bed";
    private static final String FASTA_EXT = ".fasta";

    private static final String FLD_SEQUENCE = "Sequence";
    private static final String FLD_QUALITY_SCORE = "QualityScore";
    private static final String FLD_GC_CONTENT = "GCContent";
    private static final String FLD_TARGET_TYPE = "TargetType";
    private static final String FLD_TARGET_EXTRA_INFO = "TargetExtra";
    private static final String FLD_REJECT_REASON = "RejectReason";
    private static final String FLD_TARGET_START = "TargetPositionStart";
    private static final String FLD_TARGET_END = "TargetPositionEnd";
    private static final String FLD_EVAL_CRITERIA = "EvalCriteria";

    private static final List<String> PANEL_PROBES_COLUMNS = List.of(
            FLD_CHROMOSOME, FLD_POSITION_START, FLD_POSITION_END, FLD_SEQUENCE,
            FLD_TARGET_TYPE, FLD_TARGET_EXTRA_INFO,
            FLD_QUALITY_SCORE, FLD_GC_CONTENT);

    private static final List<String> REJECTED_REGIONS_COLUMNS = List.of(
            FLD_CHROMOSOME, FLD_POSITION_START, FLD_POSITION_END,
            FLD_TARGET_TYPE, FLD_TARGET_EXTRA_INFO,
            FLD_REJECT_REASON);

    private static final List<String> CANDIDATE_PROBES_COLUMNS = List.of(
            FLD_CHROMOSOME, FLD_POSITION_START, FLD_POSITION_END,
            FLD_TARGET_START, FLD_TARGET_END,
            FLD_TARGET_TYPE, FLD_TARGET_EXTRA_INFO,
            FLD_QUALITY_SCORE, FLD_GC_CONTENT,
            FLD_EVAL_CRITERIA, FLD_REJECT_REASON);

    private static final int CANDIDATE_PROBES_BUFFER_SIZE = 1_000_000;

    private static final Logger LOGGER = LogManager.getLogger(OutputWriter.class);

    public OutputWriter(final String panelProbesFileStem, final String targetRegionsFilePath, final String rejectedRegionsFileStem,
            final String candidateProbesFilePath) throws IOException
    {
        mPanelProbesTsvWriter =
                new DelimFileWriter<>(panelProbesFileStem + TSV_EXT, PANEL_PROBES_COLUMNS, OutputWriter::writePanelProbesTsvRow);
        mPanelProbesBedWriter = createBufferedWriter(panelProbesFileStem + BED_EXT);
        mPanelProbesFastaWriter = createBufferedWriter(panelProbesFileStem + FASTA_EXT);

        mTargetRegionsWriter = createBufferedWriter(targetRegionsFilePath);

        mRejectedRegionsTsvWriter =
                new DelimFileWriter<>(
                        rejectedRegionsFileStem + TSV_EXT, REJECTED_REGIONS_COLUMNS, OutputWriter::writeRejectedRegionsTsvRow);
        mRejectedRegionsBedWriter = createBufferedWriter(rejectedRegionsFileStem + BED_EXT);

        if(candidateProbesFilePath == null)
        {
            mCandidateProbesWriter = null;
            mCandidateProbesBuffer = null;
        }
        else
        {
            mCandidateProbesWriter =
                    new DelimFileWriter<>(candidateProbesFilePath, CANDIDATE_PROBES_COLUMNS, OutputWriter::writeCandidateProbesRow);
            mCandidateProbesBuffer = new ArrayList<>(CANDIDATE_PROBES_BUFFER_SIZE);
        }
    }

    public void writePanelProbes(List<Probe> probes) throws IOException
    {
        LOGGER.debug("Writing {} panel probes to file", probes.size());

        // Must be sorted for BED files since some tools expect sorted order.
        probes = probes.stream().sorted(Comparator.comparing(Probe::region)).toList();

        for(Probe probe : probes)
        {
            if(!probe.accepted())
            {
                // If this happens there's a code bug.
                throw new IllegalArgumentException("Shouldn't be writing rejected probes");
            }

            mPanelProbesTsvWriter.writeRow(probe);
            writePanelProbesBedRow(probe);
            writePanelProbesFastaRecord(probe);
        }
    }

    private static void writePanelProbesTsvRow(final Probe probe, DelimFileWriter.Row row)
    {
        row.set(FLD_CHROMOSOME, probe.region().chromosome());
        row.set(FLD_POSITION_START, probe.region().start());
        row.set(FLD_POSITION_END, probe.region().end());
        row.set(FLD_SEQUENCE, probe.sequence());
        row.set(FLD_QUALITY_SCORE, probe.qualityScore());
        row.set(FLD_GC_CONTENT, probe.gcContent());
        row.set(FLD_TARGET_TYPE, probe.target().metadata().type().name());
        row.set(FLD_TARGET_EXTRA_INFO, probe.target().metadata().extra());
    }

    private void writePanelProbesBedRow(final Probe probe) throws IOException
    {
        mPanelProbesBedWriter.write(formatBedRow(probe.region(), probeBedName(probe)));
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
        mTargetRegionsWriter.write(formatBedRow(region.region(), targetMetadataToBedName(region.metadata())));
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
        row.set(FLD_TARGET_TYPE, region.target().metadata().type().name());
        row.set(FLD_TARGET_EXTRA_INFO, region.target().metadata().extra());
        row.set(FLD_REJECT_REASON, region.reason());
    }

    private void writeRejectedRegionsBedRow(final RejectedRegion region) throws IOException
    {
        mRejectedRegionsBedWriter.write(formatBedRow(region.region(), targetMetadataToBedName(region.target().metadata())));
    }

    public void writeCandidateProbe(final Probe probe)
    {
        if(mCandidateProbesWriter != null)
        {
            // Buffer probes to improve performance.
            mCandidateProbesBuffer.add(probe);
            checkFlushCandidateProbes(false);
        }
    }

    private void checkFlushCandidateProbes(boolean force)
    {
        if(mCandidateProbesBuffer.size() >= CANDIDATE_PROBES_BUFFER_SIZE || force)
        {
            writeCandidateProbes(mCandidateProbesBuffer);
            mCandidateProbesBuffer.clear();
        }
    }

    private void writeCandidateProbes(final List<Probe> probes)
    {
        LOGGER.debug("Writing {} candidate probes to file", probes.size());
        probes.forEach(mCandidateProbesWriter::writeRow);
    }

    private static void writeCandidateProbesRow(final Probe probe, DelimFileWriter.Row row)
    {
        // To save disk space, don't write these fields:
        //   - Probe end position (implied from start and probe length)
        //   - Base sequence (implied from probe region)
        row.set(FLD_CHROMOSOME, probe.region().chromosome());
        row.set(FLD_POSITION_START, probe.region().start());
        row.set(FLD_TARGET_START, probe.target().region().start());
        row.set(FLD_TARGET_END, probe.target().region().end());
        row.set(FLD_TARGET_TYPE, probe.target().metadata().type().name());
        row.set(FLD_TARGET_EXTRA_INFO, probe.target().metadata().extra());
        row.setOrNull(FLD_QUALITY_SCORE, probe.qualityScore());
        row.setOrNull(FLD_GC_CONTENT, probe.gcContent());
        row.setOrNull(FLD_EVAL_CRITERIA, probe.evalCriteria().toString());
        row.setOrNull(FLD_REJECT_REASON, probe.rejectionReason());
    }

    private static String probeBedName(final Probe probe)
    {
        // Purposely unbox here to throw on nulls.
        double qualityScore = probe.qualityScore();
        double gcContent = probe.gcContent();
        String baseName = targetMetadataToBedName(probe.target().metadata());
        return format("%s:QS=%.2f:GC=%.2f", baseName, qualityScore, gcContent);
    }

    private static String targetMetadataToBedName(final TargetMetadata info)
    {
        return format("%s:%s", info.type().name(), info.extra());
    }

    private static String formatBedRow(final ChrBaseRegion region, String name)
    {
        return format("%s\t%d\t%d\t%s\n", region.chromosome(), region.start() - 1, region.end(), name);
    }

    private static String getProbeLabel(final Probe probe)
    {
        // This should be unique enough.
        TargetMetadata metadata = probe.target().metadata();
        return format("%s:%s:%d", metadata.type().name(), metadata.extra(), probe.region().start());
    }

    @Override
    public void close() throws IOException
    {
        LOGGER.debug("Flushing and closing output files");

        mPanelProbesTsvWriter.close();
        mPanelProbesBedWriter.close();
        mPanelProbesFastaWriter.close();

        mTargetRegionsWriter.close();

        mRejectedRegionsTsvWriter.close();
        mRejectedRegionsBedWriter.close();

        if(mCandidateProbesWriter != null)
        {
            checkFlushCandidateProbes(true);
            mCandidateProbesWriter.close();
        }
    }
}
