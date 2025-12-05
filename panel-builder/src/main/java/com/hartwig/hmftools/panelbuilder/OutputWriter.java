package com.hartwig.hmftools.panelbuilder;

import static java.lang.String.format;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BED_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CANDIDATE_PROBES_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CANDIDATE_TARGET_REGIONS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.FASTA_EXTENSION;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_STATS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PANEL_PROBES_FILE_STEM;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_TARGETED_REGIONS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.REJECTED_FEATURES_FILE_STEM;
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
    private final DelimFileWriter<RejectedFeature> mRejectedFeaturesTsvWriter;
    private final BufferedWriter mRejectedFeaturesBedWriter;
    private final BufferedWriter mCandidateTargetRegionsBedWriter;
    @Nullable
    private final DelimFileWriter<Probe> mCandidateProbesTsvWriter;
    @Nullable
    private final ArrayList<Probe> mCandidateProbesBuffer;
    private final DelimFileWriter<Genes.GeneStats> mGeneStatsTsvWriter;
    private final DelimFileWriter<SampleVariants.VariantInfo> mSampleVariantInfoTsvWriter;
    private int mProbeId = 0;

    private enum PanelProbesColumns
    {
        StartRegion,
        InsertSequence,
        EndRegion,
        Sequence,
        TargetedStart,
        TargetedEnd,
        TargetType,
        TargetExtra,
        QualityScore,
        GCContent
    }

    private enum RejectedFeaturesColumns
    {
        Region,
        ProbeSequence,
        ProbeQualityScore,
        ProbeGCContent,
        TargetType,
        TargetExtra
    }

    private enum CandidateProbesColumns
    {
        StartRegion,
        InsertSequence,
        EndRegion,
        Sequence,
        TargetType,
        TargetExtra,
        QualityScore,
        GCContent,
        EvalCriteria
    }

    private enum GeneStatsColumns
    {
        GeneName,
        ProbeCount
    }

    private enum SampleVariantInfoColumns
    {
        Variant,
        FilterReason
    }

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

        String panelProbesTsvFile = outputFilePath.apply(PANEL_PROBES_FILE_STEM + TSV_EXTENSION);
        String panelProbesBedFile = outputFilePath.apply(PANEL_PROBES_FILE_STEM + BED_EXTENSION);
        String panelProbesFastaFile = outputFilePath.apply(PANEL_PROBES_FILE_STEM + FASTA_EXTENSION);
        String probeTargetedRegionsBedFile = outputFilePath.apply(PROBE_TARGETED_REGIONS_FILE_NAME);
        String rejectedFeaturesTsvFile = outputFilePath.apply(REJECTED_FEATURES_FILE_STEM + TSV_EXTENSION);
        String rejectedFeaturesBedFile = outputFilePath.apply(REJECTED_FEATURES_FILE_STEM + BED_EXTENSION);
        String candidateTargetRegionsBedFile = outputFilePath.apply(CANDIDATE_TARGET_REGIONS_FILE_NAME);
        String candidateProbesTsvFile = outputFilePath.apply(CANDIDATE_PROBES_FILE_NAME);
        String geneStatsTsvFile = outputFilePath.apply(GENE_STATS_FILE_NAME);
        String sampleVariantInfoTsvFile = outputFilePath.apply(SAMPLE_VARIANT_INFO_FILE_NAME);

        mPanelProbesTsvWriter =
                new DelimFileWriter<>(panelProbesTsvFile, PanelProbesColumns.values(), OutputWriter::writePanelProbesTsvRow);
        mPanelProbesBedWriter = createBufferedWriter(panelProbesBedFile);
        mPanelProbesFastaWriter = createBufferedWriter(panelProbesFastaFile);

        mProbeTargetedRegionsBedWriter = createBufferedWriter(probeTargetedRegionsBedFile);

        mRejectedFeaturesTsvWriter =
                new DelimFileWriter<>(rejectedFeaturesTsvFile, RejectedFeaturesColumns.values(), OutputWriter::writeRejectedFeaturesTsvRow);
        mRejectedFeaturesBedWriter = createBufferedWriter(rejectedFeaturesBedFile);

        mCandidateTargetRegionsBedWriter = createBufferedWriter(candidateTargetRegionsBedFile);

        if(verboseOutput)
        {
            mCandidateProbesTsvWriter =
                    new DelimFileWriter<>(candidateProbesTsvFile, CandidateProbesColumns.values(), OutputWriter::writeCandidateProbesRow);
            mCandidateProbesBuffer = new ArrayList<>(CANDIDATE_PROBES_BUFFER_SIZE);
        }
        else
        {
            mCandidateProbesTsvWriter = null;
            mCandidateProbesBuffer = null;
        }

        mGeneStatsTsvWriter = new DelimFileWriter<>(geneStatsTsvFile, GeneStatsColumns.values(), OutputWriter::writeGeneStatsRow);

        mSampleVariantInfoTsvWriter =
                new DelimFileWriter<>(sampleVariantInfoTsvFile, SampleVariantInfoColumns.values(), OutputWriter::writeSampleVariantInfoRow);
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
        // TODO? write region orientation too
        SequenceDefinition definition = probe.definition();
        ChrBaseRegion start = definition.startRegion();
        ChrBaseRegion end = definition.endRegion();
        row.setOrNull(PanelProbesColumns.StartRegion, start == null ? null : start.toString());
        row.set(PanelProbesColumns.InsertSequence, definition.insertSequence());
        row.setOrNull(PanelProbesColumns.EndRegion, end == null ? null : end.toString());
        row.setOrNull(PanelProbesColumns.Sequence, probe.sequence());
        row.set(PanelProbesColumns.TargetedStart, probe.targetedRange().startOffset());
        row.set(PanelProbesColumns.TargetedEnd, probe.targetedRange().endOffset());
        row.setOrNull(PanelProbesColumns.QualityScore, probe.qualityScore());
        row.setOrNull(PanelProbesColumns.GCContent, probe.gcContent());
        row.set(PanelProbesColumns.TargetType, probe.metadata().type().name());
        row.set(PanelProbesColumns.TargetExtra, probe.metadata().extraInfo());
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

    public void writeRejectedFeatures(List<RejectedFeature> rejectedFeatures) throws IOException
    {
        LOGGER.debug("Writing {} rejected features to file", rejectedFeatures.size());

        // Must be sorted for BED files since some tools expect sorted order.
        rejectedFeatures = rejectedFeatures.stream()
                .sorted(Comparator.comparing(RejectedFeature::region, Comparator.nullsLast(Comparator.naturalOrder()))).toList();

        for(RejectedFeature rejectedFeature : rejectedFeatures)
        {
            mRejectedFeaturesTsvWriter.writeRow(rejectedFeature);
            if(rejectedFeature.region() != null)
            {
                writeRejectedFeaturesBedRow(rejectedFeature);
            }
        }
    }

    private static void writeRejectedFeaturesTsvRow(final RejectedFeature rejectedFeature, DelimFileWriter.Row row)
    {
        row.setOrNull(RejectedFeaturesColumns.Region, rejectedFeature.region() == null ? null : rejectedFeature.region().toString());
        row.setOrNull(RejectedFeaturesColumns.ProbeSequence, rejectedFeature.probe() == null ? null : rejectedFeature.probe().sequence());
        row.setOrNull(RejectedFeaturesColumns.ProbeQualityScore, rejectedFeature.probe() == null ? null : rejectedFeature.probe().qualityScore());
        row.setOrNull(RejectedFeaturesColumns.ProbeGCContent, rejectedFeature.probe() == null ? null : rejectedFeature.probe().gcContent());
        row.set(RejectedFeaturesColumns.TargetType, rejectedFeature.metadata().type().name());
        row.set(RejectedFeaturesColumns.TargetExtra, rejectedFeature.metadata().extraInfo());
    }

    private void writeRejectedFeaturesBedRow(final RejectedFeature rejectedFeature) throws IOException
    {
        mRejectedFeaturesBedWriter.write(
                formatBedRow(requireNonNull(rejectedFeature.region()), targetMetadataToBedName(rejectedFeature.metadata())));
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
        row.setOrNull(CandidateProbesColumns.StartRegion, start == null ? null : start.toString());
        row.set(CandidateProbesColumns.InsertSequence, definition.insertSequence());
        row.setOrNull(CandidateProbesColumns.EndRegion, end == null ? null : end.toString());
        row.setOrNull(CandidateProbesColumns.Sequence, probe.sequence());
        row.set(CandidateProbesColumns.TargetType, probe.metadata().type().name());
        row.set(CandidateProbesColumns.TargetExtra, probe.metadata().extraInfo());
        row.setOrNull(CandidateProbesColumns.QualityScore, probe.qualityScore());
        row.setOrNull(CandidateProbesColumns.GCContent, probe.gcContent());
        row.set(CandidateProbesColumns.EvalCriteria, requireNonNull(probe.evaluationCriteria()).toString());
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
        row.set(GeneStatsColumns.GeneName, stats.geneName());
        row.set(GeneStatsColumns.ProbeCount, stats.probeCount());
    }

    public void writeSampleVariantInfos(final List<SampleVariants.VariantInfo> variantInfos)
    {
        variantInfos.forEach(mSampleVariantInfoTsvWriter::writeRow);
    }

    private static void writeSampleVariantInfoRow(final SampleVariants.VariantInfo variantInfo, DelimFileWriter.Row row)
    {
        row.set(SampleVariantInfoColumns.Variant, variantInfo.variant());
        row.set(SampleVariantInfoColumns.FilterReason, variantInfo.filterReason() == null ? "PASS" : variantInfo.filterReason());
    }

    @Override
    public void close() throws IOException
    {
        LOGGER.debug("Flushing and closing output files");

        mPanelProbesTsvWriter.close();
        mPanelProbesBedWriter.close();
        mPanelProbesFastaWriter.close();

        mProbeTargetedRegionsBedWriter.close();

        mRejectedFeaturesTsvWriter.close();
        mRejectedFeaturesBedWriter.close();

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
