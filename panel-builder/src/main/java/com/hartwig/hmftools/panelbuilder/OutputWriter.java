package com.hartwig.hmftools.panelbuilder;

import static java.lang.String.format;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CANDIDATE_PROBES_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CANDIDATE_TARGET_REGIONS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_STATS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PANEL_PROBES_FILE_STEM;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_TARGETED_REGIONS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.REJECTED_REGIONS_FILE_STEM;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_VARIANT_INFO_FILE_NAME;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.function.Function;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.panelbuilder.samplevariants.SampleVariants;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

// Writes all file output data.
public class OutputWriter implements AutoCloseable
{
    private final DelimFileWriter<Probe> mPanelProbesTsvWriter;
    private final BufferedWriter mPanelProbesBedWriter;
    private final BufferedWriter mPanelProbesFastaWriter;
    private final BufferedWriter mProbeTargetedRegionsBedWriter;
    private final DelimFileWriter<RejectedRegion> mRejectedRegionsTsvWriter;
    private final BufferedWriter mRejectedRegionsBedWriter;
    private final BufferedWriter mCandidateTargetRegionsBedWriter;
    @Nullable
    private final DelimFileWriter<Probe> mCandidateProbesTsvWriter;
    @Nullable
    private final ArrayList<Probe> mCandidateProbesBuffer;
    private final DelimFileWriter<Genes.GeneStats> mGeneStatsTsvWriter;
    private final DelimFileWriter<SampleVariants.VariantInfo> mSampleVariantInfoTsvWriter;
    private int mProbeId = 0;

    private static final String TSV_EXT = ".tsv";
    private static final String BED_EXT = ".bed";
    private static final String FASTA_EXT = ".fasta";

    // TODO: switch to enum columns
    private static final String FLD_START_REGION = "StartRegion";
    private static final String FLD_INSERT_SEQ = "InsertSequence";
    private static final String FLD_END_REGION = "EndRegion";
    private static final String FLD_SEQUENCE = "Sequence";
    private static final String FLD_TARGETED_START = "TargetedStart";
    private static final String FLD_TARGETED_END = "TargetedEnd";
    private static final String FLD_QUALITY_SCORE = "QualityScore";
    private static final String FLD_GC_CONTENT = "GCContent";
    private static final String FLD_TARGET_TYPE = "TargetType";
    private static final String FLD_TARGET_EXTRA_INFO = "TargetExtra";
    private static final String FLD_EVAL_CRITERIA = "EvalCriteria";
    private static final String FLD_PROBE_COUNT = "ProbeCount";
    private static final String FLD_VARIANT = "Variant";
    private static final String FLD_FILTER_REASON = "FilterReason";

    private static final List<String> PANEL_PROBES_COLUMNS = List.of(
            FLD_START_REGION, FLD_INSERT_SEQ, FLD_END_REGION, FLD_SEQUENCE,
            FLD_TARGETED_START, FLD_TARGETED_END,
            FLD_TARGET_TYPE, FLD_TARGET_EXTRA_INFO,
            FLD_QUALITY_SCORE, FLD_GC_CONTENT);

    private static final List<String> REJECTED_REGIONS_COLUMNS = List.of(
            FLD_CHROMOSOME, FLD_POSITION_START, FLD_POSITION_END,
            FLD_TARGET_TYPE, FLD_TARGET_EXTRA_INFO);

    private static final List<String> CANDIDATE_PROBES_COLUMNS = List.of(
            FLD_START_REGION, FLD_INSERT_SEQ, FLD_END_REGION, FLD_SEQUENCE,
            FLD_TARGET_TYPE, FLD_TARGET_EXTRA_INFO,
            FLD_QUALITY_SCORE, FLD_GC_CONTENT,
            FLD_EVAL_CRITERIA);

    private static final List<String> GENE_STATS_COLUMNS = List.of(
            FLD_GENE_NAME,
            FLD_PROBE_COUNT
    );

    private static final List<String> SAMPLE_VARIANT_INFO_COLUMNS = List.of(
            FLD_VARIANT,
            FLD_FILTER_REASON
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
        String probeTargetedRegionsBedFile = outputFilePath.apply(PROBE_TARGETED_REGIONS_FILE_NAME);
        String rejectedRegionsTsvFile = outputFilePath.apply(REJECTED_REGIONS_FILE_STEM + TSV_EXT);
        String rejectedRegionsBedFile = outputFilePath.apply(REJECTED_REGIONS_FILE_STEM + BED_EXT);
        String candidateTargetRegionsBedFile = outputFilePath.apply(CANDIDATE_TARGET_REGIONS_FILE_NAME);
        String candidateProbesTsvFile = outputFilePath.apply(CANDIDATE_PROBES_FILE_NAME);
        String geneStatsTsvFile = outputFilePath.apply(GENE_STATS_FILE_NAME);
        String sampleVariantInfoTsvFile = outputFilePath.apply(SAMPLE_VARIANT_INFO_FILE_NAME);

        mPanelProbesTsvWriter = new DelimFileWriter<>(panelProbesTsvFile, PANEL_PROBES_COLUMNS, OutputWriter::writePanelProbesTsvRow);
        mPanelProbesBedWriter = createBufferedWriter(panelProbesBedFile);
        mPanelProbesFastaWriter = createBufferedWriter(panelProbesFastaFile);

        mProbeTargetedRegionsBedWriter = createBufferedWriter(probeTargetedRegionsBedFile);

        mRejectedRegionsTsvWriter =
                new DelimFileWriter<>(rejectedRegionsTsvFile, REJECTED_REGIONS_COLUMNS, OutputWriter::writeRejectedRegionsTsvRow);
        mRejectedRegionsBedWriter = createBufferedWriter(rejectedRegionsBedFile);

        mCandidateTargetRegionsBedWriter = createBufferedWriter(candidateTargetRegionsBedFile);

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

        mSampleVariantInfoTsvWriter =
                new DelimFileWriter<>(sampleVariantInfoTsvFile, SAMPLE_VARIANT_INFO_COLUMNS, OutputWriter::writeSampleVariantInfoRow);
    }

    public void writePanelProbes(List<Probe> probes) throws IOException
    {
        LOGGER.debug("Writing {} panel probes to file", probes.size());

        // TODO: should there be a probe ID which matches between TSV, BED, and FASTA?

        // Must be sorted for BED files since some tools expect sorted order.
        probes = probes.stream().sorted(Comparator.comparing(
                probe -> probe.definition().singleRegionOrNull(), Comparator.nullsLast(Comparator.naturalOrder()))).toList();

        for(Probe probe : probes)
        {
            // A few basic checks that might reveal bugs in the code elsewhere.
            if(!probe.accepted())
            {
                throw new IllegalArgumentException("Should only be writing accepted probes");
            }
            if(probe.definition().baseLength() != PROBE_LENGTH || probe.sequence().length() != PROBE_LENGTH)
            {
                throw new IllegalArgumentException("Should only be writing probes of length " + PROBE_LENGTH);
            }

            mPanelProbesTsvWriter.writeRow(probe);
            if(probe.definition().isSingleRegion())
            {
                writePanelProbesBedRow(probe);
            }
            writePanelProbesFastaRecord(probe);
        }
    }

    private static void writePanelProbesTsvRow(final Probe probe, DelimFileWriter.Row row)
    {
        // TODO: should write region orientation too
        SequenceDefinition definition = probe.definition();
        ChrBaseRegion start = definition.startRegion();
        ChrBaseRegion end = definition.endRegion();
        row.setOrNull(FLD_START_REGION, start == null ? null : start.toString());
        row.set(FLD_INSERT_SEQ, definition.insertSequence());
        row.setOrNull(FLD_END_REGION, end == null ? null : end.toString());
        row.setOrNull(FLD_SEQUENCE, probe.sequence());
        row.set(FLD_TARGETED_START, probe.targetedRange().startOffset());
        row.set(FLD_TARGETED_END, probe.targetedRange().endOffset());
        row.setOrNull(FLD_QUALITY_SCORE, probe.qualityScore());
        row.setOrNull(FLD_GC_CONTENT, probe.gcContent());
        row.set(FLD_TARGET_TYPE, probe.metadata().type().name());
        row.set(FLD_TARGET_EXTRA_INFO, probe.metadata().extraInfo());
    }

    private void writePanelProbesBedRow(final Probe probe) throws IOException
    {
        mPanelProbesBedWriter.write(formatBedRow(probe.definition().singleRegion(), probeBedName(probe)));
    }

    private void writePanelProbesFastaRecord(final Probe probe) throws IOException
    {
        String label = probeFastaLabel(probe);
        String sequence = probe.sequence();
        if(sequence == null)
        {
            // If this happens there's a code bug.
            throw new IllegalArgumentException("Probe must have sequence data to write FASTA");
        }
        mPanelProbesFastaWriter.write(format(">%s\n%s\n", label, sequence));
    }

    public void writeProbeTargetedRegions(List<TargetRegion> regions) throws IOException
    {
        LOGGER.debug("Writing {} probe targeted regions to file", regions.size());

        // Must be sorted for BED files since some tools expect sorted order.
        regions = regions.stream().sorted(Comparator.comparing(TargetRegion::region)).toList();

        for(TargetRegion region : regions)
        {
            writeProbeTargetedRegionsBedRow(region);
        }
    }

    private void writeProbeTargetedRegionsBedRow(final TargetRegion region) throws IOException
    {
        writeTargetRegionBedRow(region, mProbeTargetedRegionsBedWriter);
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
    }

    private void writeRejectedRegionsBedRow(final RejectedRegion region) throws IOException
    {
        mRejectedRegionsBedWriter.write(formatBedRow(region.region(), targetMetadataToBedName(region.metadata())));
    }

    public void writeCandidateTargetRegions(List<TargetRegion> regions) throws IOException
    {
        LOGGER.debug("Writing {} candidate target regions to file", regions.size());

        // Must be sorted for BED files since some tools expect sorted order.
        regions = regions.stream().sorted(Comparator.comparing(TargetRegion::region)).toList();

        for(TargetRegion region : regions)
        {
            writeTargetRegionBedRow(region, mCandidateTargetRegionsBedWriter);
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
        SequenceDefinition definition = probe.definition();
        ChrBaseRegion start = definition.startRegion();
        ChrBaseRegion end = definition.endRegion();
        row.setOrNull(FLD_START_REGION, start == null ? null : start.toString());
        row.set(FLD_INSERT_SEQ, definition.insertSequence());
        row.setOrNull(FLD_END_REGION, end == null ? null : end.toString());
        row.set(FLD_SEQUENCE, probe.sequence());
        row.set(FLD_TARGET_TYPE, probe.metadata().type().name());
        row.set(FLD_TARGET_EXTRA_INFO, probe.metadata().extraInfo());
        row.setOrNull(FLD_QUALITY_SCORE, probe.qualityScore());
        row.setOrNull(FLD_GC_CONTENT, probe.gcContent());
        row.setOrNull(FLD_EVAL_CRITERIA, requireNonNull(probe.evaluationCriteria()).toString());
    }

    private static void writeTargetRegionBedRow(final TargetRegion region, BufferedWriter writer) throws IOException
    {
        writer.write(formatBedRow(region.region(), targetMetadataToBedName(region.metadata())));
    }

    private static String probeBedName(final Probe probe)
    {
        double qualityScore = requireNonNull(probe.qualityScore());
        double gcContent = requireNonNull(probe.gcContent());
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

    private String probeFastaLabel(final Probe probe)
    {
        int id = nextProbeId();
        TargetMetadata metadata = probe.metadata();
        return format("probe%d:%s:%s", id, metadata.type().name(), metadata.extraInfo());
    }

    private int nextProbeId()
    {
        int id = mProbeId;
        ++mProbeId;
        return id;
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

    public void writeSampleVariantInfos(final List<SampleVariants.VariantInfo> variantInfos)
    {
        variantInfos.forEach(mSampleVariantInfoTsvWriter::writeRow);
    }

    private static void writeSampleVariantInfoRow(final SampleVariants.VariantInfo variantInfo, DelimFileWriter.Row row)
    {
        row.set(FLD_VARIANT, variantInfo.variant());
        row.set(FLD_FILTER_REASON, variantInfo.filterReason());
    }

    @Override
    public void close() throws IOException
    {
        LOGGER.debug("Flushing and closing output files");

        mPanelProbesTsvWriter.close();
        mPanelProbesBedWriter.close();
        mPanelProbesFastaWriter.close();

        mProbeTargetedRegionsBedWriter.close();

        mRejectedRegionsTsvWriter.close();
        mRejectedRegionsBedWriter.close();

        mCandidateTargetRegionsBedWriter.close();

        if(mCandidateProbesTsvWriter != null)
        {
            checkFlushCandidateProbes(true);
            mCandidateProbesTsvWriter.close();
        }

        mGeneStatsTsvWriter.close();

        mSampleVariantInfoTsvWriter.close();
    }
}
