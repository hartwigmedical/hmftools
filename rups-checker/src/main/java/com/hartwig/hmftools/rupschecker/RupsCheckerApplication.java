package com.hartwig.hmftools.rupschecker;

import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.slicing.SlicerFactory;
import com.hartwig.hmftools.common.variant.SomaticTruthSetVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.predicate.VariantFilter;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.common.variant.vcf.VCFSomaticFile;
import com.hartwig.hmftools.common.variant.vcf.VCFSomaticTruthSetFile;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import javax.xml.stream.XMLStreamException;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class RupsCheckerApplication {

    private static final Logger LOGGER = LogManager.getLogger(RupsCheckerApplication.class);

    private static final String SOMATIC_MELTED_VCF = "somatic_melted_vcf";
    private static final String GERMLINE_MELTED_VCF = "germline_melted_vcf";
    private static final String TRUTH_SET_BED = "truth_set_bed";
    private static final String TRUTH_SET_VCF = "truth_set_vcf";
    private static final String PIPELINE_RUN = "run_dir";
    private static final String MD_OUTPUT = "md_output";
    private static final String SOMATIC = "somatic";
    private static final String GERMLINE = "germline";

    public static void main(final String... args)
            throws ParseException, IOException, XMLStreamException, HartwigException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        //UnrecognizedOptionException 
        final String somaticMeltedVCF = cmd.getOptionValue(SOMATIC_MELTED_VCF);
        final String truthSetBED = cmd.getOptionValue(TRUTH_SET_BED);
        final String truthSetVCF = cmd.getOptionValue(TRUTH_SET_VCF);
        final String pipelineRun = cmd.getOptionValue(PIPELINE_RUN);

        boolean produceSomaticStatistics = cmd.hasOption(SOMATIC);
        boolean produceGermlineStatistics = cmd.hasOption(GERMLINE);
        if (!(produceSomaticStatistics || produceGermlineStatistics)) {
            produceSomaticStatistics = true;
            produceGermlineStatistics = true;
        }

        boolean showHelpInfo = false;

        if (produceSomaticStatistics) {
            if (null == somaticMeltedVCF && null == pipelineRun) {
                showHelpInfo = true;
                LOGGER.error("No Somatic variant data file nor a pipeline run specified.");
            }
            if (null == truthSetBED && null == pipelineRun) {
                showHelpInfo = true;
                LOGGER.error("No Truth set BED file or a pipeline run specified.");
            }
            if (null == truthSetVCF && null == pipelineRun) {
                showHelpInfo = true;
                LOGGER.error("No Truth set VCF file or a pipeline run specified.");
            }
        }

        final String germlineMeltedVCF = cmd.getOptionValue(GERMLINE_MELTED_VCF);
        if (produceGermlineStatistics) {
            if (null == germlineMeltedVCF && null == pipelineRun) {
                showHelpInfo = true;
                LOGGER.error("No Germline variant data file nor a pipeline run specified.");
            }
        }

        if (showHelpInfo) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Run Precesion and Sensitivity Checker (RUPS)", options);
            System.exit(1);
        }

        new RupsCheckerApplication(somaticMeltedVCF,
                germlineMeltedVCF,
                truthSetBED,
                truthSetVCF,
                pipelineRun,
                produceSomaticStatistics,
                produceGermlineStatistics,
                cmd.hasOption(MD_OUTPUT)).run();
    }

    @NotNull
    private final String somaticMeltedVCF;
    @NotNull
    private final String germlineMeltedVCF;
    @NotNull
    private final String truthSetBED;
    @NotNull
    private final String truthSetVCF;
    // optional
    private final String pipelineRunDir;

    private final boolean produceSomaticStatistics;
    private final boolean produceGermlineStatistics;

    private final boolean produceMDOutput;

    private final Map<GermlinePrecisionAndSensitivity, GermlinePrecisionAndSensitivity> germlineStats = new HashMap<>();
    private final Map<SomaticPrecisionAndSensitivity, SomaticPrecisionAndSensitivity> somaticStats = new HashMap<>();

    public RupsCheckerApplication(@NotNull String somaticMeltedVCF,
            @NotNull String germlineMeltedVCF,
            @NotNull String truthSetBED,
            @NotNull String truthSetVCF,
            String pipelineRunDir,
            boolean produceSomaticStatistics,
            boolean produceGermlineStatistics,
            boolean produceMDOutput) {
        this.somaticMeltedVCF = somaticMeltedVCF;
        this.germlineMeltedVCF = germlineMeltedVCF;
        this.truthSetBED = truthSetBED;
        this.truthSetVCF = truthSetVCF;
        this.pipelineRunDir = pipelineRunDir;
        this.produceSomaticStatistics = produceSomaticStatistics;
        this.produceGermlineStatistics = produceGermlineStatistics;
        this.produceMDOutput = produceMDOutput;
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();

        options.addOption(SOMATIC_MELTED_VCF, true, "Somatic melted VCF file.");
        options.addOption(GERMLINE_MELTED_VCF, true, "Germline melted VCF file.");
        options.addOption(TRUTH_SET_BED, true, "Truth Set BED file.");
        options.addOption(TRUTH_SET_VCF, true, "Truth Set VCF file.");
        options.addOption(PIPELINE_RUN, true, "Directory holding a HMF Somatic pipeline run.");
        options.addOption(MD_OUTPUT, false, "Produces table ouput in MD format, ASCII is default.");
        options.addOption(SOMATIC, false, "Produce somatic statistics. Default is somatic and germline.");
        options.addOption(GERMLINE, false, "Produce germline statistics. Default is somatic and germline.");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private void run() throws IOException, HartwigException {
        LOGGER.info("RUPS - RUn Precision and Sensitivity Checker");
        LOGGER.info("--------------------------------------------");

        runSomaticStatistics();
        runGermlineStatistics();
        displayStatistics();

        LOGGER.info("RUPS - Done");
    }

    private void runSomaticStatistics() throws IOException, HartwigException {
        if (produceSomaticStatistics) {
            LOGGER.info("  -> Loading somatic variant data: " + somaticMeltedVCF);
            final VCFSomaticFile somaticVariantFile = VCFFileLoader.loadSomaticVCF(somaticMeltedVCF);
            LOGGER.info("  -> Filtering somatic variant data for PASS only");
            final List<SomaticVariant> somaticVariants = VariantFilter.passOnly(somaticVariantFile.variants());  
            if (somaticVariants.isEmpty()) {
                LOGGER.warn("  -> No somatic variant data. Nothing to do.");
                return;
            }
            
            LOGGER.info("  -> Loading truthSet BED data: " + truthSetBED);
            final Slicer truthSetSlicer = SlicerFactory.fromBedFile(truthSetBED);
            LOGGER.info("  -> Filtering somatic variants on TruthSet bed file");
//        List<SomaticVariant> truthSetFilteredVariants = new ArrayList<>();    
//        somaticVariants.stream().filter((variant) -> (truthSetSlicer.includes(variant))).forEachOrdered((variant) -> {
//            truthSetFilteredVariants.add(variant);
//        });
            LOGGER.info("  -> Loading truthSet VCF data: " + truthSetVCF);
            final VCFSomaticTruthSetFile truthSetFile = VCFFileLoader.loadSomaticTruthSetVCF(truthSetVCF);
            LOGGER.info("  -> Filtering truthSet for PASS only");
            final List<SomaticTruthSetVariant> truthSetVariants = VariantFilter.passOnly(truthSetFile.variants());
            
            class SomaticVariantIndex {
                private final String chromosome;
                private final long position;

                public SomaticVariantIndex(String chromosome, long position) {
                    this.chromosome = chromosome;
                    this.position = position;
                }

                @Override
                public int hashCode() {
                    int hash = 5;
                    hash = 71 * hash + Objects.hashCode(this.chromosome);
                    hash = 71 * hash + (int) (this.position ^ (this.position >>> 32));
                    return hash;
                }

                @Override
                public boolean equals(Object obj) {
                    if (this == obj) {
                        return true;
                    }
                    if (obj == null) {
                        return false;
                    }
                    if (getClass() != obj.getClass()) {
                        return false;
                    }
                    final SomaticVariantIndex other = (SomaticVariantIndex) obj;
                    if (this.position != other.position) {
                        return false;
                    }
                    if (!Objects.equals(this.chromosome, other.chromosome)) {
                        return false;
                    }
                    return true;
                }
            }
            
            Map<SomaticVariantIndex, SomaticTruthSetVariant> somaticTruthSetVariantsLookup = new HashMap<>();
            long truthSetSNP = 0;
            long truthSetINDEL = 0;
            long truthSetUNDEF = 0;
            for (SomaticTruthSetVariant truthSetSomaticVariant : truthSetVariants) {
                switch (truthSetSomaticVariant.type()) {
                    case INDEL:
                        truthSetINDEL++;
                        break;
                    case SNP:
                        truthSetSNP++;
                        break;
                    default:
                        truthSetUNDEF++;
                }
                somaticTruthSetVariantsLookup.put(new SomaticVariantIndex(truthSetSomaticVariant.chromosome(), truthSetSomaticVariant.position()), truthSetSomaticVariant);
            }
            LOGGER.info("  -> truthSetSNP: " + truthSetSNP);
            LOGGER.info("  -> truthSetINDEL: " + truthSetINDEL);
            LOGGER.warn("  -> truthSetUNDEF: " + truthSetUNDEF);
            
            LOGGER.info("  -> Checking somatic variants agains truth set");
            final Set<String> allCallers = new HashSet<>();
            Map<SomaticVariantIndex, SomaticVariant> somaticVariantsLookup = new HashMap<>();
            for (SomaticVariant somaticVariant : somaticVariants) {
                final SomaticVariantIndex key = new SomaticVariantIndex(somaticVariant.chromosome(), somaticVariant.position());
                if (somaticTruthSetVariantsLookup.containsKey(key)) {
                    //SomaticTruthSetVariant truthSetVariant = somaticTruthSetVariantsLookup.get(key);
                    for (String caller: somaticVariant.callers()) {
                        final SomaticPrecisionAndSensitivity rupsKey = new SomaticPrecisionAndSensitivity(somaticVariant.type(), caller);
                        if (null == somaticStats.get(rupsKey)) {
                            somaticStats.put(rupsKey, rupsKey);
                            allCallers.add(caller);
                        }
                        rupsKey.incTruePositive();
                    }
                } else {
                    // Not in the truthset, for all callers we have a false positive
                    for (String caller : somaticVariant.callers()) {
                        final SomaticPrecisionAndSensitivity rupsKey = new SomaticPrecisionAndSensitivity(somaticVariant.type(), caller);
                        if (null == somaticStats.get(rupsKey)) {
                            somaticStats.put(rupsKey, rupsKey);
                            allCallers.add(caller);
                        }
                        rupsKey.incFalsePositives();
                    }
                }
                somaticVariantsLookup.put(new SomaticVariantIndex(somaticVariant.chromosome(), somaticVariant.position()), somaticVariant);
            }
            for (SomaticTruthSetVariant somaticTruthSetVariant : truthSetVariants) {
                final SomaticVariantIndex key = new SomaticVariantIndex(somaticTruthSetVariant.chromosome(), somaticTruthSetVariant.position());
                if (!somaticVariantsLookup.containsKey(key)) {
                    // in the truth set but not called by any of our callers
                    for (String caller : allCallers) {
                        final SomaticPrecisionAndSensitivity rupsKey = new SomaticPrecisionAndSensitivity(somaticTruthSetVariant.type(), caller);
                        if (null == somaticStats.get(rupsKey)) {
                            somaticStats.put(rupsKey, rupsKey);
                        }
                        rupsKey.incFalseNegatives();
                    }
                }
            }
        }
    }

    private void runGermlineStatistics() throws IOException, HartwigException {
        if (produceGermlineStatistics) {
            LOGGER.info("  -> Loading germline data from: " + germlineMeltedVCF);
            // TODO. Howto read a Germline truth set?? 
            //        final VCFGermlineFile germlineVariantFile = VCFFileLoader.loadGermlineVCF(germlineMeltedVCF);
            //        final List<GermlineVariant> germlineVariants = VariantFilter.passOnly(germlineVariantFile.variants());
        }
    }
    
       private void displayStatistics() {
        if (produceGermlineStatistics) {
            LOGGER.info("");
            String headerMarker = produceMDOutput ? "**" : "";
            LOGGER.info(headerMarker + "Germline precision & sensitivity" + headerMarker);
            LOGGER.info("");
            boolean headerShown = false;
            for (GermlinePrecisionAndSensitivity germlineStat : germlineStats.keySet()) {
                if (!headerShown) {
                    headerShown = true;
                    LOGGER.info(germlineStat.getHeader(produceMDOutput));
                }
                LOGGER.info(germlineStat.getValueString(produceMDOutput));
            }
        }

        if (produceSomaticStatistics) {
            LOGGER.info("");
            String headerMarker = produceMDOutput ? "**" : "";
            LOGGER.info(headerMarker + "Somatic precision & sentivity" + headerMarker);
            LOGGER.info("");
            boolean headerShown = false;
            for (SomaticPrecisionAndSensitivity somaticStat : somaticStats.keySet()) {
                if (!headerShown) {
                    headerShown = true;
                    LOGGER.info(somaticStat.getHeader(produceMDOutput));
                }
                LOGGER.info(somaticStat.getValueString(produceMDOutput));
            }
        }
    }
}
