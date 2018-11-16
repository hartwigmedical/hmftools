package com.hartwig.hmftools.hotspotcaller;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.hotspot.HotspotEvidence;
import com.hartwig.hmftools.common.hotspot.HotspotEvidenceType;
import com.hartwig.hmftools.common.hotspot.ImmutableHotspotEvidence;
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

public class HotspotCallerApplication2 {

    private static final Logger LOGGER = LogManager.getLogger(HotspotCallerApplication2.class);

    private static final String HOTSPOT = "hotspot";
    private static final String TUMOR_BAM = "tumor_bam";
    private static final String REFERENCE_BAM = "reference_bam";
    private static final String INFRAME_BED = "inframe_bed";
    private static final String REF_GENOME = "ref_genome";

    private static final int PHRED_OFFSET = 33;

    public static void main(String[] args) throws IOException, ParseException {

        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String hotspotPath = cmd.getOptionValue(HOTSPOT);
        final String tumorBam = cmd.getOptionValue(TUMOR_BAM);
        final String referenceBam = cmd.getOptionValue(REFERENCE_BAM);
        final String refGenome = cmd.getOptionValue(REF_GENOME);
        final String inframeBed = cmd.getOptionValue(INFRAME_BED);

        if (hotspotPath == null || tumorBam == null || referenceBam == null || inframeBed == null || refGenome == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("HotspotCallerApplication2", options);
            System.exit(1);
        }

        LOGGER.info("Loading inframe bed regions");
        final Collection<GenomeRegion> codingRegions = BEDFileLoader.fromBedFile(inframeBed).values();

        final SamReader tumorReader = SamReaderFactory.makeDefault().open(new File(tumorBam));
        final SamReader referenceReader = SamReaderFactory.makeDefault().open(new File(referenceBam));

        LOGGER.info("Starting");
        final IndexedFastaSequenceFile refSequence = new IndexedFastaSequenceFile(new File(refGenome));

        LOGGER.info("Loading known hotspots from {}", hotspotPath);
        final Set<VariantHotspot> knownHotspots = Sets.newHashSet(VariantHotspotFile.read(hotspotPath).values());

        LOGGER.info("Looking for potential inframe indel locations ");
        final Set<VariantHotspot> inframeIndelsHotspots = new InframeIndelHotspots(codingRegions, refSequence).findInframeIndels(tumorReader);

        final Set<VariantHotspot> allHotspots = Sets.newHashSet();
        allHotspots.addAll(knownHotspots);
        allHotspots.addAll(inframeIndelsHotspots);

        LOGGER.info("Looking for evidence of hotspots in tumor bam {}", tumorBam);
        final VariantHotspotEvidenceFactory tumorEvidenceFactory = new VariantHotspotEvidenceFactory(codingRegions, refSequence, tumorReader);
        final List<VariantHotspotEvidence> tumorEvidence = tumorEvidenceFactory.evidence(allHotspots);

        LOGGER.info("Looking for evidence of hotspots in reference bam {}", tumorBam);
        final VariantHotspotEvidenceFactory referenceEvidenceFactory = new VariantHotspotEvidenceFactory(codingRegions, refSequence, referenceReader);
        final List<VariantHotspotEvidence> referenceEvidence = tumorEvidenceFactory.evidence(allHotspots);


        final List<HotspotEvidence> evidence = Lists.newArrayList();


        //        for (final VariantHotspot variantHotspot : knownHotspots.values()) {
//            LOGGER.info("LOOKING FOR HOT {}", variantHotspot);
//            final VariantHotspotEvidence tumorEvidence = tumorEvidenceFactory.evidence(variantHotspot);
//            if (tumorEvidence.altSupport() > 0) {
//                final VariantHotspotEvidence referenceEvidence = referenceEvidenceFactory.evidence(variantHotspot);
//                evidence.add(createEvidence(tumorEvidence, referenceEvidence));
//            }
//        }
//        LOGGER.info("Found {} knownHotspots in tumor", evidence.size());



        LOGGER.info("Complete");
    }

    private static HotspotEvidence createEvidence(@NotNull final VariantHotspotEvidence tumor,
            @NotNull final VariantHotspotEvidence normal) {
        return ImmutableHotspotEvidence.builder()
                .from(tumor)
                .ref(tumor.ref())
                .alt(tumor.alt())
                .qualityScore(tumor.altQuality())
                .tumorAltCount(tumor.altSupport())
                .tumorRefCount(tumor.refSupport())
                .tumorReads(tumor.readDepth())
                .normalRefCount(normal.refSupport())
                .normalAltCount(normal.altQuality())
                .normalReads(normal.readDepth())
                .type(HotspotEvidenceType.KNOWN)
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
        options.addOption(TUMOR_BAM, true, "Tumor bam file.");
        options.addOption(REFERENCE_BAM, true, "Reference bam file.");
        options.addOption(HOTSPOT, true, "Hotspot input file.");
        options.addOption(INFRAME_BED, true, "Hotspot input file.");
        options.addOption(REF_GENOME, true, "Hotspot input file.");
        return options;
    }
}
