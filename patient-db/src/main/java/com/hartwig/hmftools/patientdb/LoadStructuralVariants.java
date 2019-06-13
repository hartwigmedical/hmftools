package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.getValueNotNull;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFile;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFileLoader;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class LoadStructuralVariants {

    private static final Logger LOGGER = LogManager.getLogger(LoadStructuralVariants.class);

    private static final String TUMOR_SAMPLE = "tumor";
    private static final String VCF = "structural_vcf";
    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";
    private static final String REF_GENOME = "ref_genome";
    private static final String ALIAS = "alias";
    private static final String SV_DATA_DIRECTORY = "sv_data_dir";

    public static void main(@NotNull final String[] args) throws ParseException, IOException, SQLException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final DatabaseAccess dbAccess = databaseAccess(cmd);

        final String tumorSample = cmd.getOptionValue(TUMOR_SAMPLE);
        final String vcfPath = cmd.getOptionValue(VCF);
        final String svDataOutputDir = cmd.getOptionValue(SV_DATA_DIRECTORY);

        LOGGER.info("Reading data from {}", vcfPath);
        final IndexedFastaSequenceFile sequenceFile =
                cmd.hasOption(REF_GENOME) ? new IndexedFastaSequenceFile(new File(cmd.getOptionValue(REF_GENOME))) : null;

        final EnrichedStructuralVariantFactory factory = new EnrichedStructuralVariantFactory(sequenceFile);
        final List<StructuralVariant> variants = StructuralVariantFileLoader.fromFile(vcfPath, new AlwaysPassFilter());

        LOGGER.info("Enriching variants");
        final List<EnrichedStructuralVariant> enrichedVariants = factory.enrich(variants);

        // generate a unique ID for each SV record
        int svId = 0;

        List<StructuralVariantData> svDataList = Lists.newArrayList();

        for (EnrichedStructuralVariant var : enrichedVariants) {
            svDataList.add(convertSvData(var, svId++));
        }

        LOGGER.info("Persisting {} SVs to db", svDataList.size());
        dbAccess.writeStructuralVariants(cmd.getOptionValue(ALIAS, tumorSample), svDataList);

        if (svDataOutputDir != null) {
            // write data to file
            try {
                final String svFilename = StructuralVariantFile.generateFilename(svDataOutputDir, tumorSample);
                StructuralVariantFile.write(svFilename, svDataList);
            } catch (IOException e) {
                LOGGER.error("failed to write SV data: {}", e.toString());
            }
        }

        LOGGER.info("Complete");
    }

    public static StructuralVariantData convertSvData(final EnrichedStructuralVariant var, int svId) {
        return ImmutableStructuralVariantData.builder()
                .id(svId)
                .startChromosome(var.chromosome(true))
                .endChromosome(var.end() == null ? "0" : var.chromosome(false))
                .startPosition(var.position(true))
                .endPosition(var.end() == null ? -1 : var.position(false))
                .startOrientation(var.orientation(true))
                .endOrientation(var.end() == null ? (byte) 0 : var.orientation(false))
                .startHomologySequence(var.start().homology())
                .endHomologySequence(var.end() == null ? "" : var.end().homology())
                .ploidy(getValueNotNull(var.ploidy()))
                .startAF(getValueNotNull(var.start().alleleFrequency()))
                .endAF(var.end() == null ? 0 : getValueNotNull(var.end().alleleFrequency()))
                .adjustedStartAF(getValueNotNull(var.start().adjustedAlleleFrequency()))
                .adjustedEndAF(var.end() == null ? 0 : getValueNotNull(var.end().adjustedAlleleFrequency()))
                .adjustedStartCopyNumber(getValueNotNull(var.start().adjustedCopyNumber()))
                .adjustedEndCopyNumber(var.end() == null ? 0 : getValueNotNull(var.end().adjustedCopyNumber()))
                .adjustedStartCopyNumberChange(getValueNotNull(var.start().adjustedCopyNumberChange()))
                .adjustedEndCopyNumberChange(var.end() == null ? 0 : getValueNotNull(var.end().adjustedCopyNumberChange()))
                .insertSequence(var.insertSequence())
                .type(var.type())
                .filter(var.filter())
                .imprecise(var.imprecise())
                .qualityScore(getValueNotNull(var.qualityScore()))
                .event(getValueNotNull(var.event()))
                .startTumorVariantFragmentCount(getValueNotNull(var.start().tumorVariantFragmentCount()))
                .startTumorReferenceFragmentCount(getValueNotNull(var.start().tumorReferenceFragmentCount()))
                .startNormalVariantFragmentCount(getValueNotNull(var.start().normalVariantFragmentCount()))
                .startNormalReferenceFragmentCount(getValueNotNull(var.start().normalReferenceFragmentCount()))
                .endTumorVariantFragmentCount(var.end() == null ? 0 : getValueNotNull(var.end().tumorVariantFragmentCount()))
                .endTumorReferenceFragmentCount(var.end() == null ? 0 : getValueNotNull(var.end().tumorReferenceFragmentCount()))
                .endNormalVariantFragmentCount(var.end() == null ? 0 : getValueNotNull(var.end().normalVariantFragmentCount()))
                .endNormalReferenceFragmentCount(var.end() == null ? 0 : getValueNotNull(var.end().normalReferenceFragmentCount()))
                .startIntervalOffsetStart(getValueNotNull(var.start().startOffset()))
                .startIntervalOffsetEnd(getValueNotNull(var.start().endOffset()))
                .endIntervalOffsetStart(var.end() == null ? 0 : getValueNotNull(var.end().startOffset()))
                .endIntervalOffsetEnd(var.end() == null ? 0 : getValueNotNull(var.end().endOffset()))
                .inexactHomologyOffsetStart(getValueNotNull(var.start().inexactHomologyOffsetStart()))
                .inexactHomologyOffsetEnd(getValueNotNull(var.start().inexactHomologyOffsetEnd()))
                .startLinkedBy(getValueNotNull(var.startLinkedBy()))
                .endLinkedBy(getValueNotNull(var.endLinkedBy()))
                .vcfId(getValueNotNull(var.id()))
                .startRefContext(getValueNotNull(var.start().refGenomeContext()))
                .endRefContext(var.end() == null ? "" : getValueNotNull(var.end().refGenomeContext()))
                .recovered(var.recovered())
                .recoveryMethod((getValueNotNull(var.recoveryMethod())))
                .recoveryFilter(getValueNotNull(var.recoveryFilter()))
                .insertSequenceAlignments(getValueNotNull(var.insertSequenceAlignments()))
                .insertSequenceRepeatClass(getValueNotNull(var.insertSequenceRepeatClass()))
                .insertSequenceRepeatType(getValueNotNull(var.insertSequenceRepeatType()))
                .insertSequenceRepeatOrientation(getValueNotNull(var.insertSequenceRepeatOrientation()))
                .insertSequenceRepeatCoverage(getValueNotNull(var.insertSequenceRepeatCoverage()))
                .startAnchoringSupportDistance(var.start().anchoringSupportDistance())
                .endAnchoringSupportDistance(var.end() == null ? 0 : var.end().anchoringSupportDistance())
                .build();
    }

    @NotNull
    private static Options createBasicOptions() {
        final Options options = new Options();
        options.addOption(TUMOR_SAMPLE, true, "Name of the tumor sample. This should correspond to the value used in PURPLE.");
        options.addOption(VCF, true, "Path to the PURPLE structural variant VCF file.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(REF_GENOME, true, "Path to the (indexed) ref genome fasta file.");
        options.addOption(ALIAS, true, "Overwrite the sample name with specified alias when writing to db");
        options.addOption(SV_DATA_DIRECTORY, true, "Optional: directory to write SV data in TSV format");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static DatabaseAccess databaseAccess(@NotNull final CommandLine cmd) throws SQLException {
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);  //e.g. mysql://localhost:port/database";
        final String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }
}
