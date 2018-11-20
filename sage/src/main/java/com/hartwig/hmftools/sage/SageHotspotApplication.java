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
import com.hartwig.hmftools.common.hotspot.HotspotEvidence;
import com.hartwig.hmftools.common.hotspot.HotspotEvidenceType;
import com.hartwig.hmftools.common.hotspot.HotspotEvidenceVCF;
import com.hartwig.hmftools.common.hotspot.ImmutableHotspotEvidence;
import com.hartwig.hmftools.common.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.hotspot.InframeIndelHotspots;
import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.hotspot.VariantHotspotEvidence;
import com.hartwig.hmftools.common.hotspot.VariantHotspotEvidenceFactory;
import com.hartwig.hmftools.common.hotspot.VariantHotspotFile;
import com.hartwig.hmftools.common.region.BEDFileLoader;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SageHotspotApplication {

    private static final Logger LOGGER = LogManager.getLogger(SageHotspotApplication.class);

    private static final String OUT_PATH = "out";
    private static final String HOTSPOT = "hotspot";
    private static final String TUMOR = "tumor";
    private static final String TUMOR_BAM = "tumor_bam";
    private static final String REFERENCE = "reference";
    private static final String REFERENCE_BAM = "reference_bam";
    private static final String INFRAME_BED = "inframe_bed";
    private static final String REF_GENOME = "ref_genome";

    public static void main(String[] args) throws IOException, ParseException {

        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String hotspotPath = cmd.getOptionValue(HOTSPOT);
        final String tumorBam = cmd.getOptionValue(TUMOR_BAM);
        final String referenceBam = cmd.getOptionValue(REFERENCE_BAM);
        final String refGenome = cmd.getOptionValue(REF_GENOME);
        final String inframeBed = cmd.getOptionValue(INFRAME_BED);
        final String outputVCF = cmd.getOptionValue(OUT_PATH);
        final String referenceSample = cmd.getOptionValue(REFERENCE);
        final String tumorSample = cmd.getOptionValue(TUMOR);

        if (hotspotPath == null || tumorBam == null || referenceBam == null || inframeBed == null || refGenome == null || outputVCF == null
                || tumorSample == null || referenceSample == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SAGE", options);
            System.exit(1);
        }

        final IndexedFastaSequenceFile refSequence = new IndexedFastaSequenceFile(new File(refGenome));
        final SamReader tumorReader = SamReaderFactory.makeDefault().open(new File(tumorBam));
        final SamReader referenceReader = SamReaderFactory.makeDefault().open(new File(referenceBam));

        LOGGER.info("Loading inframe bed regions from {}", inframeBed);
        final Collection<GenomeRegion> codingRegions = BEDFileLoader.fromBedFile(inframeBed).values();

        LOGGER.info("Loading known hotspots from {}", hotspotPath);
        final Set<VariantHotspot> knownHotspots = Sets.newHashSet(VariantHotspotFile.read(hotspotPath).values());

        LOGGER.info("Looking for potential inframe indel locations ");
        final Set<VariantHotspot> allHotspots = Sets.newHashSet();
        allHotspots.addAll(knownHotspots);
        allHotspots.addAll(new InframeIndelHotspots(codingRegions, refSequence).findInframeIndels(tumorReader));

        LOGGER.info("Looking for evidence of hotspots in tumor bam {}", tumorBam);
        final VariantHotspotEvidenceFactory tumorEvidenceFactory =
                new VariantHotspotEvidenceFactory(codingRegions, refSequence, tumorReader);
        final Map<VariantHotspot, VariantHotspotEvidence> tumorEvidence = asMap(tumorEvidenceFactory.evidence(allHotspots));

        LOGGER.info("Looking for evidence of hotspots in reference bam {}", tumorBam);
        final VariantHotspotEvidenceFactory referenceEvidenceFactory =
                new VariantHotspotEvidenceFactory(codingRegions, refSequence, referenceReader);
        final Map<VariantHotspot, VariantHotspotEvidence> referenceEvidence = asMap(referenceEvidenceFactory.evidence(allHotspots));

        final List<HotspotEvidence> evidence = Lists.newArrayList();
        for (Map.Entry<VariantHotspot, VariantHotspotEvidence> entry : tumorEvidence.entrySet()) {
            final VariantHotspot variant = entry.getKey();
            final VariantHotspotEvidence tumor = entry.getValue();
            final VariantHotspotEvidence normal = referenceEvidence.get(variant);
            if (tumor.altSupport() > 0) {
                evidence.add(createEvidence(knownHotspots.contains(variant), tumor, normal));
            }
        }

        LOGGER.info("Writing output to {}", outputVCF);
        Collections.sort(evidence);
        new HotspotEvidenceVCF(referenceSample, tumorSample).write(outputVCF, evidence);

        LOGGER.info("Complete");
    }

    @NotNull
    private static Map<VariantHotspot, VariantHotspotEvidence> asMap(@NotNull final List<VariantHotspotEvidence> evidence) {
        return evidence.stream().collect(Collectors.toMap(x -> ImmutableVariantHotspotImpl.builder().from(x).build(), x -> x));
    }

    private static HotspotEvidence createEvidence(boolean known, @NotNull final VariantHotspotEvidence tumor,
            @NotNull final VariantHotspotEvidence normal) {
        return ImmutableHotspotEvidence.builder()
                .from(tumor)
                .ref(tumor.ref())
                .alt(tumor.alt())
                .qualityScore(tumor.altQuality())
                .tumorAltCount(tumor.altSupport())
                .tumorRefCount(tumor.refSupport())
                .tumorReads(tumor.readDepth())
                .normalAltCount(normal.altSupport())
                .normalRefCount(normal.refSupport())
                .normalReads(normal.readDepth())
                .normalIndelCount(normal.indelSupport())
                .type(known ? HotspotEvidenceType.KNOWN : HotspotEvidenceType.INFRAME)
                .build();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(OUT_PATH, true, "Tumor bam file.");
        options.addOption(TUMOR_BAM, true, "Tumor bam file.");
        options.addOption(REFERENCE_BAM, true, "Reference bam file.");
        options.addOption(HOTSPOT, true, "Hotspot input file.");
        options.addOption(INFRAME_BED, true, "Hotspot input file.");
        options.addOption(REF_GENOME, true, "Hotspot input file.");
        options.addOption(TUMOR, true, "Tumor sample name.");
        options.addOption(REFERENCE, true, "Reference sample name.");
        return options;
    }

}
