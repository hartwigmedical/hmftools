package com.hartwig.hmftools.panelbuilder;

import static java.lang.String.format;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BED_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.FASTA_EXTENSION;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.mergeOverlapAndAdjacentRegions;
import static com.hartwig.hmftools.panelbuilder.Utils.combineStringUnique;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Comparator;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.BiConsumer;

// Writes the output files for one probe panel (DNA or RNA): probes (TSV/FASTA/BED), covered target regions, covered regions, candidate
// target regions, rejected features, and gene stats. DNA and RNA share this logic; they differ only in the probes TSV layout (DNA uses the
// start/end BasicProbeLayout, limited to two regions; RNA lists all regions generically to allow spliced probes) and the FASTA id prefix.
public class ProbeOutputWriter implements AutoCloseable
{
    private final boolean mRna;
    private final DelimFileWriter<Probe> mProbesTsvWriter;
    private final BufferedWriter mProbesBedWriter;
    private final BufferedWriter mProbesFastaWriter;
    private final BufferedWriter mCoveredTargetRegionsBedWriter;
    private final BufferedWriter mCoveredRegionsBedWriter;
    private final DelimFileWriter<RejectedFeature> mRejectedFeaturesTsvWriter;
    private final BufferedWriter mRejectedFeaturesBedWriter;
    private final BufferedWriter mCandidateTargetRegionsBedWriter;
    private final DelimFileWriter<GeneStats> mGeneStatsTsvWriter;
    private int mProbeId = 0;

    private enum ProbesColumns
    {
        StartRegion,
        StartRegionOrient,
        MiddleSequence,
        EndRegion,
        EndRegionOrient,
        Sequence,
        TargetedStart,
        TargetedEnd,
        TargetType,
        TargetExtra,
        QualityScore,
        GCContent
    }

    // RNA probes may span more than two regions (spliced across exon junctions), so unlike the DNA probe format all segments are listed
    // generically. Ref segments are genome-forward (see RNA design decision).
    private enum RnaProbesColumns
    {
        Segments,
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
        TargetExtra,
        RejectionReason,
    }

    private enum GeneStatsColumns
    {
        GeneName,
        ProbeCount
    }

    private static final Logger LOGGER = LogManager.getLogger(ProbeOutputWriter.class);

    public ProbeOutputWriter(final Function<String, String> outputFilePath, final String probesStem, final String coveredTargetRegionsFile,
            final String coveredRegionsFile, final String rejectedFeaturesStem, final String candidateTargetRegionsFile,
            final String geneStatsFile, boolean rna) throws IOException
    {
        mRna = rna;

        Enum<?>[] probesColumns = rna ? RnaProbesColumns.values() : ProbesColumns.values();
        BiConsumer<Probe, DelimFileWriter.Row> probesRowWriter = rna ? this::writeRnaProbesRow : this::writeProbesRow;
        mProbesTsvWriter = new DelimFileWriter<>(outputFilePath.apply(probesStem + TSV_EXTENSION), probesColumns, probesRowWriter);
        mProbesBedWriter = createBufferedWriter(outputFilePath.apply(probesStem + BED_EXTENSION));
        mProbesFastaWriter = createBufferedWriter(outputFilePath.apply(probesStem + FASTA_EXTENSION));

        mCoveredTargetRegionsBedWriter = createBufferedWriter(outputFilePath.apply(coveredTargetRegionsFile));
        mCoveredRegionsBedWriter = createBufferedWriter(outputFilePath.apply(coveredRegionsFile));

        mRejectedFeaturesTsvWriter = new DelimFileWriter<>(
                outputFilePath.apply(rejectedFeaturesStem + TSV_EXTENSION), RejectedFeaturesColumns.values(),
                ProbeOutputWriter::writeRejectedFeaturesTsvRow);
        mRejectedFeaturesBedWriter = createBufferedWriter(outputFilePath.apply(rejectedFeaturesStem + BED_EXTENSION));

        mCandidateTargetRegionsBedWriter = createBufferedWriter(outputFilePath.apply(candidateTargetRegionsFile));

        mGeneStatsTsvWriter =
                new DelimFileWriter<>(outputFilePath.apply(geneStatsFile), GeneStatsColumns.values(), ProbeOutputWriter::writeGeneStatsRow);
    }

    public void writeProbes(final List<Probe> probes) throws IOException
    {
        LOGGER.debug("Writing {} panel probes to file", probes.size());

        // Sort for consistent output. RNA uses the full SequenceDefinition (a total order over multi-region probes); DNA sorts by its single
        // region (variant probes with no single region sort last).
        Comparator<Probe> order = mRna
                ? Comparator.comparing(Probe::definition)
                : Comparator.comparing(probe -> probe.definition().singleRegionOrNull(), Comparator.nullsLast(Comparator.naturalOrder()));
        List<Probe> probesSorted = probes.stream().sorted(order).toList();

        for(Probe probe : probesSorted)
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

            mProbesTsvWriter.writeRow(probe);
            writeFastaRecord(probe);
        }

        // Must be sorted for BED files since some tools expect sorted order.
        List<NamedRegion> bedRegions = probes.stream()
                .flatMap(probe ->
                {
                    String name = probeBedName(probe);
                    return probe.definition().regions().stream().map(region -> new NamedRegion(region, name));
                })
                .sorted()
                .toList();
        for(NamedRegion region : bedRegions)
        {
            writeBedRow(region, mProbesBedWriter);
        }
    }

    private void writeProbesRow(final Probe probe, DelimFileWriter.Row row)
    {
        // Note this throws if the probe is not in the right format.
        // The DNA panel probes output format represents at most a start and end region. Multi-region (e.g. spliced) probes require a
        // different output; RNA has its own format.
        BasicProbeLayout layout = BasicProbeLayout.from(probe.definition());
        ChrBaseRegion start = layout.startRegion();
        Orientation startOrientation = layout.startOrientation();
        ChrBaseRegion end = layout.endRegion();
        Orientation endOrientation = layout.endOrientation();
        row.setOrNull(ProbesColumns.StartRegion, start == null ? null : start.toString());
        row.setOrNull(ProbesColumns.StartRegionOrient, startOrientation == null ? null : startOrientation.asChar());
        row.setOrNull(ProbesColumns.MiddleSequence, layout.insertSequence());
        row.setOrNull(ProbesColumns.EndRegion, end == null ? null : end.toString());
        row.setOrNull(ProbesColumns.EndRegionOrient, endOrientation == null ? null : endOrientation.asChar());
        row.setOrNull(ProbesColumns.Sequence, probe.sequence());
        row.set(ProbesColumns.TargetedStart, probe.targetedRange().startOffset());
        row.set(ProbesColumns.TargetedEnd, probe.targetedRange().endOffset());
        row.setOrNull(ProbesColumns.QualityScore, probe.qualityScore());
        row.setOrNull(ProbesColumns.GCContent, probe.gcContent());
        row.set(ProbesColumns.TargetType, probe.metadata().type().name());
        row.set(ProbesColumns.TargetExtra, probe.metadata().extraInfo());
    }

    private void writeRnaProbesRow(final Probe probe, DelimFileWriter.Row row)
    {
        String segments = probe.definition().segments().stream().map(ProbeOutputWriter::formatSegment).collect(Collectors.joining(";"));
        row.set(RnaProbesColumns.Segments, segments);
        row.setOrNull(RnaProbesColumns.Sequence, probe.sequence());
        row.set(RnaProbesColumns.TargetedStart, probe.targetedRange().startOffset());
        row.set(RnaProbesColumns.TargetedEnd, probe.targetedRange().endOffset());
        row.set(RnaProbesColumns.TargetType, probe.metadata().type().name());
        row.set(RnaProbesColumns.TargetExtra, probe.metadata().extraInfo());
        row.setOrNull(RnaProbesColumns.QualityScore, probe.qualityScore());
        row.setOrNull(RnaProbesColumns.GCContent, probe.gcContent());
    }

    // Serializes one probe segment: a ref segment as "chromosome:start-end:orientation", an insert segment as its literal bases.
    private static String formatSegment(final SequenceSegment segment)
    {
        if(segment instanceof RefSegment ref)
        {
            return format("%s:%s", ref.region(), ref.orientation().asChar());
        }
        else
        {
            return ((InsertSeqSegment) segment).sequence();
        }
    }

    private void writeFastaRecord(final Probe probe) throws IOException
    {
        String label = probeFastaLabel(probe);
        String sequence = probe.sequence();
        if(sequence == null)
        {
            // If this happens there's a code bug.
            throw new IllegalArgumentException("Probe must have sequence data to write FASTA");
        }
        mProbesFastaWriter.write(format(">%s\n%s\n", label, sequence));
    }

    public void writeCoveredTargetRegions(final List<TargetRegion> regions) throws IOException
    {
        writeTargetRegionsBed(regions, mCoveredTargetRegionsBedWriter);
    }

    public void writeCandidateTargetRegions(final List<TargetRegion> regions) throws IOException
    {
        writeTargetRegionsBed(regions, mCandidateTargetRegionsBedWriter);
    }

    private static void writeTargetRegionsBed(final List<TargetRegion> regions, final BufferedWriter writer) throws IOException
    {
        LOGGER.debug("Writing {} target regions to file", regions.size());

        // Must be sorted for BED files since some tools expect sorted order.
        List<TargetRegion> regionsSorted = regions.stream().sorted(Comparator.comparing(TargetRegion::region)).toList();

        for(TargetRegion region : regionsSorted)
        {
            writeTargetRegionBedRow(region, writer);
        }
    }

    public void writeCoveredRegions(final List<Probe> probes) throws IOException
    {
        List<NamedRegion> regions = createCoveredRegions(probes);

        LOGGER.debug("Writing {} covered regions to file", regions.size());

        // Must be sorted for BED files since some tools expect sorted order.
        List<NamedRegion> regionsSorted = regions.stream().sorted().toList();

        for(NamedRegion region : regionsSorted)
        {
            writeBedRow(region, mCoveredRegionsBedWriter);
        }
    }

    private static List<NamedRegion> createCoveredRegions(final List<Probe> probes)
    {
        Stream<NamedRegion> regions = probes.stream()
                .flatMap(probe -> probe.definition()
                        .regions()
                        .stream()
                        .map(region -> new NamedRegion(region, targetMetadataToBedName(probe.metadata()))));
        return mergeOverlapAndAdjacentRegions(regions, NamedRegion::region, ProbeOutputWriter::mergeCoveredRegion);
    }

    private static NamedRegion mergeCoveredRegion(final ChrBaseRegion mergedRegion, final NamedRegion r1, final NamedRegion r2)
    {
        return new NamedRegion(
                mergedRegion, combineStringUnique(r1.name(), r2.name(), (s1, s2) -> format("%s | %s", r1.name(), r2.name())));
    }

    public void writeRejectedFeatures(final List<RejectedFeature> rejectedFeatures) throws IOException
    {
        LOGGER.debug("Writing {} rejected features to file", rejectedFeatures.size());

        // Must be sorted for BED files since some tools expect sorted order.
        List<RejectedFeature> rejectedFeaturesSorted = rejectedFeatures.stream()
                .sorted(Comparator.comparing(RejectedFeature::region, Comparator.nullsLast(Comparator.naturalOrder()))).toList();

        for(RejectedFeature rejectedFeature : rejectedFeaturesSorted)
        {
            mRejectedFeaturesTsvWriter.writeRow(rejectedFeature);
            // TODO: should write the regions from the probe too?
            if(rejectedFeature.region() != null)
            {
                writeBedRow(requireNonNull(rejectedFeature.region()), targetMetadataToBedName(rejectedFeature.metadata()),
                        mRejectedFeaturesBedWriter);
            }
        }
    }

    private static void writeRejectedFeaturesTsvRow(final RejectedFeature rejectedFeature, DelimFileWriter.Row row)
    {
        Probe probe = rejectedFeature.probe();
        row.setOrNull(RejectedFeaturesColumns.Region, rejectedFeature.region() == null ? null : rejectedFeature.region().toString());
        row.setOrNull(RejectedFeaturesColumns.ProbeSequence, probe == null ? null : probe.sequence());
        row.setOrNull(RejectedFeaturesColumns.ProbeQualityScore, probe == null ? null : probe.qualityScore());
        row.setOrNull(RejectedFeaturesColumns.ProbeGCContent, probe == null ? null : probe.gcContent());
        row.set(RejectedFeaturesColumns.TargetType, rejectedFeature.metadata().type().name());
        row.set(RejectedFeaturesColumns.TargetExtra, rejectedFeature.metadata().extraInfo());
        row.set(RejectedFeaturesColumns.RejectionReason, probe == null ? null : requireNonNull(probe.evaluationResult()).rejectionInfo());
    }

    public void writeGeneStats(final List<GeneStats> geneStats)
    {
        geneStats.forEach(mGeneStatsTsvWriter::writeRow);
    }

    private static void writeGeneStatsRow(final GeneStats stats, DelimFileWriter.Row row)
    {
        row.set(GeneStatsColumns.GeneName, stats.geneName());
        row.set(GeneStatsColumns.ProbeCount, stats.probeCount());
    }

    private static void writeTargetRegionBedRow(final TargetRegion region, BufferedWriter writer) throws IOException
    {
        writeBedRow(region.region(), targetMetadataToBedName(region.metadata()), writer);
    }

    private static void writeBedRow(final NamedRegion region, BufferedWriter writer) throws IOException
    {
        writeBedRow(region.region(), region.name(), writer);
    }

    private static void writeBedRow(final ChrBaseRegion region, final String name, BufferedWriter writer) throws IOException
    {
        writer.write(format("%s\t%d\t%d\t%s\n", region.chromosome(), region.start() - 1, region.end(), name));
    }

    private static String probeBedName(final Probe probe)
    {
        String name = targetMetadataToBedName(probe.metadata());

        Double qualityScore = probe.qualityScore();
        if(qualityScore != null)
        {
            name += format(":QS=%.2f", qualityScore);
        }

        Double gcContent = probe.gcContent();
        if(gcContent != null)
        {
            name += format(":GC=%.2f", gcContent);
        }

        return name;
    }

    private static String targetMetadataToBedName(final TargetMetadata info)
    {
        return format("%s:%s", info.type().name(), info.extraInfo());
    }

    // Probe ids are prefixed per panel (dna_/rna_) with separate counters, so ids are never ambiguous between DNA and RNA output.
    private String probeFastaLabel(final Probe probe)
    {
        int id = mProbeId++;
        TargetMetadata metadata = probe.metadata();
        return format("%s_%d:%s:%s", mRna ? "rna" : "dna", id, metadata.type().name(), metadata.extraInfo());
    }

    @Override
    public void close() throws IOException
    {
        mProbesTsvWriter.close();
        mProbesBedWriter.close();
        mProbesFastaWriter.close();
        mCoveredTargetRegionsBedWriter.close();
        mCoveredRegionsBedWriter.close();
        mRejectedFeaturesTsvWriter.close();
        mRejectedFeaturesBedWriter.close();
        mCandidateTargetRegionsBedWriter.close();
        mGeneStatsTsvWriter.close();
    }
}
