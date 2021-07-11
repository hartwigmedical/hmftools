package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;
import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.valueNotNull;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariantLeg;
import com.hartwig.hmftools.common.sv.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.DatabaseUtil;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadStructuralVariants {

    private static final Logger LOGGER = LogManager.getLogger(LoadStructuralVariants.class);

    private static final String SAMPLE = "sample";
    private static final String SV_VCF = "structural_vcf";

    public static void main(@NotNull String[] args) throws ParseException, IOException, SQLException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);
        DatabaseAccess dbAccess = databaseAccess(cmd);

        if (dbAccess == null) {
            LOGGER.error("Failed to create DB connection");
            return;
        }

        String sample = cmd.getOptionValue(SAMPLE);
        LOGGER.info("Loading structural variants for {}", sample);

        String vcfPath = cmd.getOptionValue(SV_VCF);
        LOGGER.info("Reading SV data from {}", vcfPath);
        List<StructuralVariant> variants = StructuralVariantFileLoader.fromFile(vcfPath, new AlwaysPassFilter());
        List<EnrichedStructuralVariant> enrichedVariants = new EnrichedStructuralVariantFactory().enrich(variants);

        // Generate a unique ID for each SV record
        int svId = 0;

        List<StructuralVariantData> structuralVariants = Lists.newArrayList();
        for (EnrichedStructuralVariant variant : enrichedVariants) {
            structuralVariants.add(convertSvData(variant, svId++));
        }

        LOGGER.info("Writing {} SVs", structuralVariants.size());
        dbAccess.writeStructuralVariants(sample, structuralVariants);

        LOGGER.info("Complete");
    }

    @NotNull
    private static StructuralVariantData convertSvData(@NotNull EnrichedStructuralVariant var, int svId) {
        EnrichedStructuralVariantLeg start = var.start();
        EnrichedStructuralVariantLeg end = var.end();
        String filter = var.filter();
        Boolean imprecise = var.imprecise();

        return ImmutableStructuralVariantData.builder()
                .id(svId)
                .startChromosome(start.chromosome())
                .endChromosome(end == null ? "0" : end.chromosome())
                .startPosition((int)start.position())
                .endPosition(end == null ? -1 : (int)end.position())
                .startOrientation(start.orientation())
                .endOrientation(end == null ? (byte) 0 : end.orientation())
                .startHomologySequence(start.homology())
                .endHomologySequence(end == null ? "" : end.homology())
                .junctionCopyNumber(DatabaseUtil.valueNotNull(var.junctionCopyNumber()))
                .startAF(DatabaseUtil.valueNotNull(start.alleleFrequency()))
                .endAF(end == null ? 0 : DatabaseUtil.valueNotNull(end.alleleFrequency()))
                .adjustedStartAF(DatabaseUtil.valueNotNull(start.adjustedAlleleFrequency()))
                .adjustedEndAF(end == null ? 0 : DatabaseUtil.valueNotNull(end.adjustedAlleleFrequency()))
                .adjustedStartCopyNumber(DatabaseUtil.valueNotNull(start.adjustedCopyNumber()))
                .adjustedEndCopyNumber(end == null ? 0 : DatabaseUtil.valueNotNull(end.adjustedCopyNumber()))
                .adjustedStartCopyNumberChange(DatabaseUtil.valueNotNull(start.adjustedCopyNumberChange()))
                .adjustedEndCopyNumberChange(end == null ? 0 : DatabaseUtil.valueNotNull(end.adjustedCopyNumberChange()))
                .insertSequence(var.insertSequence())
                .type(var.type())
                .filter(filter != null ? filter : "Unknown")
                .imprecise(imprecise != null ? imprecise : false)
                .qualityScore(DatabaseUtil.valueNotNull(var.qualityScore()))
                .event(valueNotNull(var.event()))
                .startTumorVariantFragmentCount(DatabaseUtil.valueNotNull(start.tumorVariantFragmentCount()))
                .startTumorReferenceFragmentCount(DatabaseUtil.valueNotNull(start.tumorReferenceFragmentCount()))
                .startNormalVariantFragmentCount(DatabaseUtil.valueNotNull(start.normalVariantFragmentCount()))
                .startNormalReferenceFragmentCount(DatabaseUtil.valueNotNull(start.normalReferenceFragmentCount()))
                .endTumorVariantFragmentCount(end == null ? 0 : DatabaseUtil.valueNotNull(end.tumorVariantFragmentCount()))
                .endTumorReferenceFragmentCount(end == null ? 0 : DatabaseUtil.valueNotNull(end.tumorReferenceFragmentCount()))
                .endNormalVariantFragmentCount(end == null ? 0 : DatabaseUtil.valueNotNull(end.normalVariantFragmentCount()))
                .endNormalReferenceFragmentCount(end == null ? 0 : DatabaseUtil.valueNotNull(end.normalReferenceFragmentCount()))
                .startIntervalOffsetStart(DatabaseUtil.valueNotNull(start.startOffset()))
                .startIntervalOffsetEnd(DatabaseUtil.valueNotNull(start.endOffset()))
                .endIntervalOffsetStart(end == null ? 0 : DatabaseUtil.valueNotNull(end.startOffset()))
                .endIntervalOffsetEnd(end == null ? 0 : DatabaseUtil.valueNotNull(end.endOffset()))
                .inexactHomologyOffsetStart(DatabaseUtil.valueNotNull(start.inexactHomologyOffsetStart()))
                .inexactHomologyOffsetEnd(DatabaseUtil.valueNotNull(start.inexactHomologyOffsetEnd()))
                .startLinkedBy(valueNotNull(var.startLinkedBy()))
                .endLinkedBy(valueNotNull(var.endLinkedBy()))
                .vcfId(valueNotNull(var.id()))
                .startRefContext(valueNotNull(start.refGenomeContext()))
                .endRefContext(end == null ? "" : valueNotNull(end.refGenomeContext()))
                .recovered(var.recovered())
                .recoveryMethod((valueNotNull(var.recoveryMethod())))
                .recoveryFilter(valueNotNull(var.recoveryFilter()))
                .insertSequenceAlignments(valueNotNull(var.insertSequenceAlignments()))
                .insertSequenceRepeatClass(valueNotNull(var.insertSequenceRepeatClass()))
                .insertSequenceRepeatType(valueNotNull(var.insertSequenceRepeatType()))
                .insertSequenceRepeatOrientation(DatabaseUtil.valueNotNull(var.insertSequenceRepeatOrientation()))
                .insertSequenceRepeatCoverage(DatabaseUtil.valueNotNull(var.insertSequenceRepeatCoverage()))
                .startAnchoringSupportDistance(start.anchoringSupportDistance())
                .endAnchoringSupportDistance(end == null ? 0 : end.anchoringSupportDistance())
                .build();
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Name of the tumor sample. This should correspond to the value used in PURPLE.");
        options.addOption(SV_VCF, true, "Path to the PURPLE structural variant VCF file.");
        addDatabaseCmdLineArgs(options);

        return options;
    }
}
