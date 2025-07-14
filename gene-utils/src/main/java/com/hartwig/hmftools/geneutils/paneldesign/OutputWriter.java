package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class OutputWriter implements AutoCloseable
{
    private final DelimFileWriter<EvaluatedProbe> mPanelProbesTsvWriter;
    private final BufferedWriter mPanelProbesBedWriter;
    private final BufferedWriter mTargetRegionsWriter;
    private final DelimFileWriter<RejectedRegion> mRejectedRegionsTsvWriter;
    private final BufferedWriter mRejectedRegionsBedWriter;
    private final DelimFileWriter<EvaluatedProbe> mCandidateProbesWriter;
    private final List<EvaluatedProbe> mCandidateProbesBuffer;

    private static final String TSV_EXT = ".tsv";
    private static final String BED_EXT = ".bed";

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

    public OutputWriter(final String panelProbesFile, final String targetRegionsFile, final String rejectedRegionsFile,
            final String candidateProbesFile) throws IOException
    {
        mPanelProbesTsvWriter =
                new DelimFileWriter<>(panelProbesFile + TSV_EXT, PANEL_PROBES_COLUMNS, OutputWriter::writePanelProbesTsvRow);
        mPanelProbesBedWriter = createBufferedWriter(panelProbesFile + BED_EXT);
        mTargetRegionsWriter = createBufferedWriter(targetRegionsFile);
        mRejectedRegionsTsvWriter =
                new DelimFileWriter<>(rejectedRegionsFile + TSV_EXT, REJECTED_REGIONS_COLUMNS, OutputWriter::writeRejectedRegionsTsvRow);
        mRejectedRegionsBedWriter = createBufferedWriter(rejectedRegionsFile + BED_EXT);
        if(candidateProbesFile == null)
        {
            mCandidateProbesWriter = null;
            mCandidateProbesBuffer = null;
        }
        else
        {
            mCandidateProbesWriter =
                    new DelimFileWriter<>(candidateProbesFile, CANDIDATE_PROBES_COLUMNS, OutputWriter::writeCandidateProbesRow);
            mCandidateProbesBuffer = new ArrayList<>(CANDIDATE_PROBES_BUFFER_SIZE);
        }
    }

    public void writePanelProbes(final List<EvaluatedProbe> probes) throws IOException
    {
        LOGGER.debug("Writing {} panel probes to file", probes.size());
        for(EvaluatedProbe probe : probes)
        {
            mPanelProbesTsvWriter.writeRow(probe);
            writePanelProbesBedRow(probe);
        }
    }

    private static void writePanelProbesTsvRow(final EvaluatedProbe probe, DelimFileWriter.Row row)
    {
        row.set(FLD_CHROMOSOME, probe.candidate().probeRegion().chromosome());
        row.set(FLD_POSITION_START, probe.candidate().probeRegion().start());
        row.set(FLD_POSITION_END, probe.candidate().probeRegion().end());
        row.set(FLD_SEQUENCE, probe.sequence());
        row.set(FLD_QUALITY_SCORE, probe.qualityScore());
        row.set(FLD_GC_CONTENT, probe.gcContent());
        row.set(FLD_TARGET_TYPE, probe.candidate().target().metadata().type().name());
        row.set(FLD_TARGET_EXTRA_INFO, probe.candidate().target().metadata().extra());
    }

    private void writePanelProbesBedRow(final EvaluatedProbe probe) throws IOException
    {
        CandidateProbe candidate = probe.candidate();
        mPanelProbesBedWriter.write(formatBedRow(candidate.probeRegion(), targetMetadataToBedName(candidate.target().metadata())));
    }

    public void writeTargetRegions(final List<TargetRegion> regions) throws IOException
    {
        LOGGER.debug("Writing {} target regions to file", regions.size());
        for(TargetRegion region : regions)
        {
            writeTargetRegionsBedRow(region);
        }
    }

    private void writeTargetRegionsBedRow(final TargetRegion region) throws IOException
    {
        mTargetRegionsWriter.write(formatBedRow(region.region(), targetMetadataToBedName(region.metadata())));
    }

    public void writeRejectedRegions(final List<RejectedRegion> regions) throws IOException
    {
        LOGGER.debug("Writing {} rejected regions to file", regions.size());
        for(RejectedRegion region : regions)
        {
            mRejectedRegionsTsvWriter.writeRow(region);
            writeRejectedRegionsBedRow(region);
        }
    }

    private static void writeRejectedRegionsTsvRow(final RejectedRegion region, DelimFileWriter.Row row)
    {
        row.set(FLD_CHROMOSOME, region.baseRegion().chromosome());
        row.set(FLD_POSITION_START, region.baseRegion().start());
        row.set(FLD_POSITION_END, region.baseRegion().end());
        row.set(FLD_TARGET_TYPE, region.target().metadata().type().name());
        row.set(FLD_TARGET_EXTRA_INFO, region.target().metadata().extra());
        row.set(FLD_REJECT_REASON, region.reason());
    }

    private void writeRejectedRegionsBedRow(final RejectedRegion region) throws IOException
    {
        mRejectedRegionsBedWriter.write(formatBedRow(region.baseRegion(), targetMetadataToBedName(region.target().metadata())));
    }

    public void writeCandidateProbe(final EvaluatedProbe probe)
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

    private void writeCandidateProbes(final List<EvaluatedProbe> probes)
    {
        LOGGER.debug("Writing {} candidate probes to file", probes.size());
        probes.forEach(mCandidateProbesWriter::writeRow);
    }

    private static void writeCandidateProbesRow(final EvaluatedProbe probe, DelimFileWriter.Row row)
    {
        // To save disk space, don't write these fields:
        //   - Probe end position (implied from start and probe length)
        //   - Base sequence (implied from probe region)
        row.set(FLD_CHROMOSOME, probe.candidate().probeRegion().chromosome());
        row.set(FLD_POSITION_START, probe.candidate().probeRegion().start());
        row.set(FLD_TARGET_START, probe.candidate().target().region().start());
        row.set(FLD_TARGET_END, probe.candidate().target().region().end());
        row.set(FLD_TARGET_TYPE, probe.candidate().target().metadata().type().name());
        row.set(FLD_TARGET_EXTRA_INFO, probe.candidate().target().metadata().extra());
        row.setOrNull(FLD_QUALITY_SCORE, probe.qualityScore());
        row.setOrNull(FLD_GC_CONTENT, probe.gcContent());
        row.setOrNull(FLD_EVAL_CRITERIA, probe.criteria().toString());
        row.setOrNull(FLD_REJECT_REASON, probe.rejectionReason());
    }

    private static String targetMetadataToBedName(final TargetMetadata info)
    {
        return format("%s:%s", info.type().name(), info.extra());
    }

    private static String formatBedRow(final ChrBaseRegion region, String name)
    {
        return format("%s\t%d\t%d\t%s\n", region.chromosome(), region.start() - 1, region.end(), name);
    }

    @Override
    public void close()
    {
        mPanelProbesTsvWriter.close();
        mRejectedRegionsTsvWriter.close();
        if(mCandidateProbesWriter != null)
        {
            checkFlushCandidateProbes(true);
            mCandidateProbesWriter.close();
        }
    }
}
