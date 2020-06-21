package com.hartwig.hmftools.sage;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.BEDFileLoader;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.HotspotEvidence;
import com.hartwig.hmftools.common.variant.hotspot.HotspotEvidenceType;
import com.hartwig.hmftools.common.variant.hotspot.HotspotEvidenceVCF;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableHotspotEvidence;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.InframeIndelHotspots;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotEvidence;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotEvidenceFactory;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SageHotspotApplication implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(SageHotspotApplication.class);

    public static void main(String[] args) throws IOException {
        final Options options = SageHotspotApplicationConfig.createOptions();
        try (final SageHotspotApplication application = new SageHotspotApplication(options, args)) {
            application.run();
        } catch (ParseException e) {
            LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SageHotspotApplication", options);
            System.exit(1);
        }
    }

    private final SamReader tumorReader;
    private final SamReader referenceReader;
    private final SageHotspotApplicationConfig config;

    private SageHotspotApplication(final Options options, final String... args) throws ParseException {
        final CommandLine cmd = createCommandLine(args, options);
        config = SageHotspotApplicationConfig.createConfig(cmd);

        final File refGenomeFile = new File(config.refGenomePath());
        tumorReader = SamReaderFactory.makeDefault().referenceSequence(refGenomeFile).open(new File(config.tumorBamPath()));
        referenceReader = SamReaderFactory.makeDefault().referenceSequence(refGenomeFile).open(new File(config.referenceBamPath()));
    }

    private void run() throws IOException {
        final String hotspotPath = config.knownHotspotPath();
        final String tumorBam = config.tumorBamPath();
        final String referenceBam = config.referenceBamPath();
        final String codingRegionBedFile = config.codingRegionBedPath();
        final String outputVCF = config.outputFile();
        final String referenceSample = config.normal();
        final String tumorSample = config.tumor();
        final int minMappingQuality = config.minMappingQuality();
        final int minBaseQuality = config.minBaseQuality();

        final File refGenomeFile = new File(config.refGenomePath());
        final IndexedFastaSequenceFile refSequence = new IndexedFastaSequenceFile(refGenomeFile);

        LOGGER.info("Loading coding regions from {}", codingRegionBedFile);
        final Collection<GenomeRegion> codingRegions = BEDFileLoader.fromBedFile(codingRegionBedFile).values();

        LOGGER.info("Loading known hotspots from {}", hotspotPath);
        final Set<VariantHotspot> knownHotspots = Sets.newHashSet(VariantHotspotFile.read(hotspotPath).values());

        LOGGER.info("Looking for potential inframe indel locations ");
        final Set<VariantHotspot> allHotspots = Sets.newHashSet();
        allHotspots.addAll(knownHotspots);
        allHotspots.addAll(new InframeIndelHotspots(minMappingQuality, codingRegions, refSequence).findInframeIndels(tumorReader));

        LOGGER.info("Looking for evidence of hotspots in tumor bam {}", tumorBam);
        final VariantHotspotEvidenceFactory hotspotEvidenceFactory =
                new VariantHotspotEvidenceFactory(minMappingQuality, minBaseQuality, allHotspots);
        final Map<VariantHotspot, VariantHotspotEvidence> tumorEvidence = asMap(hotspotEvidenceFactory.evidence(refSequence, tumorReader));

        LOGGER.info("Looking for evidence of hotspots in reference bam {}", referenceBam);
        final Map<VariantHotspot, VariantHotspotEvidence> referenceEvidence =
                asMap(hotspotEvidenceFactory.evidence(refSequence, referenceReader));

        final List<HotspotEvidence> evidence = Lists.newArrayList();
        for (Map.Entry<VariantHotspot, VariantHotspotEvidence> entry : tumorEvidence.entrySet()) {
            final VariantHotspot variant = entry.getKey();
            final VariantHotspotEvidence tumor = entry.getValue();
            final VariantHotspotEvidence normal = referenceEvidence.get(variant);
            evidence.add(createEvidence(knownHotspots.contains(variant), tumor, normal));
        }

        LOGGER.info("Writing output to {}", outputVCF);
        Collections.sort(evidence);
        new HotspotEvidenceVCF(referenceSample,
                tumorSample,
                config.maxHetBinomialLikelihood(),
                config.minTumorReads(),
                config.minSnvVAF(),
                config.minIndelVAF(),
                config.minSnvQuality(),
                config.minIndelQuality()).write(outputVCF, evidence);
    }

    @NotNull
    private static Map<VariantHotspot, VariantHotspotEvidence> asMap(@NotNull final List<VariantHotspotEvidence> evidence) {
        return evidence.stream().collect(Collectors.toMap(x -> ImmutableVariantHotspotImpl.builder().from(x).build(), x -> x));
    }

    private static HotspotEvidence createEvidence(boolean known, @NotNull final VariantHotspotEvidence tumor,
            @Nullable final VariantHotspotEvidence normal) {
        return ImmutableHotspotEvidence.builder()
                .from(tumor)
                .ref(tumor.ref())
                .alt(tumor.alt())
                .qualityScore(tumor.altQuality())
                .tumorAltCount(tumor.altSupport())
                .tumorRefCount(tumor.refSupport())
                .tumorReads(tumor.readDepth())
                .normalAltCount(normal == null ? 0 : normal.altSupport())
                .normalRefCount(normal == null ? 0 : normal.refSupport())
                .normalReads(normal == null ? 0 : normal.readDepth())
                .normalIndelCount(normal == null ? 0 : normal.indelSupport())
                .type(known ? HotspotEvidenceType.KNOWN : HotspotEvidenceType.INFRAME)
                .build();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @Override
    public void close() throws IOException {
        tumorReader.close();
        referenceReader.close();
        LOGGER.info("Complete");
    }
}
