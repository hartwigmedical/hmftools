package com.hartwig.hmftools.serve.dao;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.actionability.ActionableEvents;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.immuno.ActionableHLA;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.actionability.util.ActionableFileFunctions;
import com.hartwig.hmftools.serve.cancertype.CancerTypeFactory;
import com.hartwig.hmftools.serve.extraction.KnownEvents;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodon;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.extraction.hotspot.KnownHotspot;

import static com.hartwig.hmftools.serve.database.tables.Actionablehotspot.ACTIONABLEHOTSPOT;
import static com.hartwig.hmftools.serve.database.tables.Actionablerange.ACTIONABLERANGE;
import static com.hartwig.hmftools.serve.database.tables.Actionablegene.ACTIONABLEGENE;
import static com.hartwig.hmftools.serve.database.tables.Actionablefusion.ACTIONABLEFUSION;
import static com.hartwig.hmftools.serve.database.tables.Actionablecharacteristic.ACTIONABLECHARACTERISTIC;
import static com.hartwig.hmftools.serve.database.tables.Actionablehla.ACTIONABLEHLA;
import static com.hartwig.hmftools.serve.database.tables.Knownhotspot.KNOWNHOTSPOT;
import static com.hartwig.hmftools.serve.database.tables.Knowncodon.KNOWNCODON;
import static com.hartwig.hmftools.serve.database.tables.Knownexon.KNOWNEXON;
import static com.hartwig.hmftools.serve.database.tables.Knowncopynumber.KNOWNCOPYNUMBER;
import static com.hartwig.hmftools.serve.database.tables.Knownfusionpair.KNOWNFUSIONPAIR;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep13;
import org.jooq.InsertValuesStep14;
import org.jooq.InsertValuesStep15;
import org.jooq.InsertValuesStep16;
import org.jooq.InsertValuesStep18;
import org.jooq.InsertValuesStep20;
import org.jooq.InsertValuesStep4;
import org.jooq.InsertValuesStep8;
import org.jooq.InsertValuesStep9;

public class ServeDAO {

    private static final Logger LOGGER = LogManager.getLogger(ServeDAO.class);

    @NotNull
    private final DSLContext context;

    ServeDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    public void deleteAll() {
        LOGGER.info("Deleting all data from SERVE database");
        context.deleteFrom(ACTIONABLEHOTSPOT).execute();
        context.deleteFrom(ACTIONABLERANGE).execute();
        context.deleteFrom(ACTIONABLEGENE).execute();
        context.deleteFrom(ACTIONABLEFUSION).execute();
        context.deleteFrom(ACTIONABLECHARACTERISTIC).execute();
        context.deleteFrom(ACTIONABLEHLA).execute();
        context.deleteFrom(KNOWNHOTSPOT).execute();
        context.deleteFrom(KNOWNCODON).execute();
        context.deleteFrom(KNOWNEXON).execute();
        context.deleteFrom(KNOWNFUSIONPAIR).execute();
        context.deleteFrom(KNOWNCOPYNUMBER).execute();
    }

    void write(@NotNull ActionableEvents actionableEvents, KnownEvents knownEvents) {

        deleteAll();

        Timestamp timestamp = new Timestamp(new Date().getTime());

        List<ActionableHotspot> actionableHotspots = actionableEvents.hotspots();
        for (List<ActionableHotspot> batch : Iterables.partition(actionableHotspots, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep16 inserter = context.insertInto(ACTIONABLEHOTSPOT,
                    ACTIONABLEHOTSPOT.MODIFIED,
                    ACTIONABLEHOTSPOT.CHROMOSOME,
                    ACTIONABLEHOTSPOT.POSITION,
                    ACTIONABLEHOTSPOT.REF,
                    ACTIONABLEHOTSPOT.ALT,
                    ACTIONABLEHOTSPOT.SOURCE,
                    ACTIONABLEHOTSPOT.SOURCEEVENT,
                    ACTIONABLEHOTSPOT.SOURCEURLS,
                    ACTIONABLEHOTSPOT.TREATMENT,
                    ACTIONABLEHOTSPOT.DRUGCLASSES,
                    ACTIONABLEHOTSPOT.APPLICABLECANCERTYPE,
                    ACTIONABLEHOTSPOT.APPLICABLEDOID,
                    ACTIONABLEHOTSPOT.BLACKLISTCANCERTYPES,
                    ACTIONABLEHOTSPOT.LEVEL,
                    ACTIONABLEHOTSPOT.DIRECTION,
                    ACTIONABLEHOTSPOT.EVIDENCEURLS);
            batch.forEach(entry -> addRecordHotspots(timestamp, inserter, entry));
            inserter.execute();
        }

        List<ActionableRange> actionableRanges = actionableEvents.ranges();
        for (List<ActionableRange> batch : Iterables.partition(actionableRanges, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep20 inserter = context.insertInto(ACTIONABLERANGE,
                    ACTIONABLERANGE.MODIFIED,
                    ACTIONABLERANGE.GENE,
                    ACTIONABLERANGE.TRANSCRIPT,
                    ACTIONABLERANGE.CHROMOSOME,
                    ACTIONABLERANGE.START,
                    ACTIONABLERANGE.END,
                    ACTIONABLERANGE.MUTATIONTYPE,
                    ACTIONABLERANGE.RANGETYPE,
                    ACTIONABLERANGE.RANGERANK,
                    ACTIONABLERANGE.SOURCE,
                    ACTIONABLERANGE.SOURCEEVENT,
                    ACTIONABLERANGE.SOURCEURLS,
                    ACTIONABLERANGE.TREATMENT,
                    ACTIONABLERANGE.DRUGCLASSES,
                    ACTIONABLERANGE.APPLICABLECANCERTYPE,
                    ACTIONABLERANGE.APPLICABLEDOID,
                    ACTIONABLERANGE.BLACKLISTCANCERTYPES,
                    ACTIONABLERANGE.LEVEL,
                    ACTIONABLERANGE.DIRECTION,
                    ACTIONABLERANGE.EVIDENCEURLS);
            batch.forEach(entry -> addRecordRanges(timestamp, inserter, entry));
            inserter.execute();
        }

        List<ActionableGene> actionableGenes = actionableEvents.genes();
        for (List<ActionableGene> batch : Iterables.partition(actionableGenes, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep14 inserter = context.insertInto(ACTIONABLEGENE,
                    ACTIONABLEGENE.MODIFIED,
                    ACTIONABLEGENE.GENE,
                    ACTIONABLEGENE.EVENT,
                    ACTIONABLEGENE.SOURCE,
                    ACTIONABLEGENE.SOURCEEVENT,
                    ACTIONABLEGENE.SOURCEURLS,
                    ACTIONABLEGENE.TREATMENT,
                    ACTIONABLEGENE.DRUGCLASSES,
                    ACTIONABLEGENE.APPLICABLECANCERTYPE,
                    ACTIONABLEGENE.APPLICABLEDOID,
                    ACTIONABLEGENE.BLACKLISTCANCERTYPES,
                    ACTIONABLEGENE.LEVEL,
                    ACTIONABLEGENE.DIRECTION,
                    ACTIONABLEGENE.EVIDENCEURLS);
            batch.forEach(entry -> addRecordGenes(timestamp, inserter, entry));
            inserter.execute();
        }

        List<ActionableFusion> actionableFusions = actionableEvents.fusions();
        for (List<ActionableFusion> batch : Iterables.partition(actionableFusions, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep18 inserter = context.insertInto(ACTIONABLEFUSION,
                    ACTIONABLEFUSION.MODIFIED,
                    ACTIONABLEFUSION.GENEUP,
                    ACTIONABLEFUSION.MINEXONUP,
                    ACTIONABLEFUSION.MAXEXONUP,
                    ACTIONABLEFUSION.GENEDOWN,
                    ACTIONABLEFUSION.MINEXONDOWN,
                    ACTIONABLEFUSION.MAXEXONDOWN,
                    ACTIONABLEFUSION.SOURCE,
                    ACTIONABLEFUSION.SOURCEEVENT,
                    ACTIONABLEFUSION.SOURCEURLS,
                    ACTIONABLEFUSION.TREATMENT,
                    ACTIONABLEFUSION.DRUGCLASSES,
                    ACTIONABLEFUSION.APPLICABLECANCERTYPE,
                    ACTIONABLEFUSION.APPLICABLEDOID,
                    ACTIONABLEFUSION.BLACKLISTCANCERTYPES,
                    ACTIONABLEFUSION.LEVEL,
                    ACTIONABLEFUSION.DIRECTION,
                    ACTIONABLEFUSION.EVIDENCEURLS);
            batch.forEach(entry -> addRecordFusions(timestamp, inserter, entry));
            inserter.execute();
        }

        List<ActionableCharacteristic> actionableCharacteristics = actionableEvents.characteristics();
        for (List<ActionableCharacteristic> batch : Iterables.partition(actionableCharacteristics, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep15 inserter = context.insertInto(ACTIONABLECHARACTERISTIC,
                    ACTIONABLECHARACTERISTIC.MODIFIED,
                    ACTIONABLECHARACTERISTIC.NAME,
                    ACTIONABLECHARACTERISTIC.COMPARATOR,
                    ACTIONABLECHARACTERISTIC.CUTOFF,
                    ACTIONABLECHARACTERISTIC.SOURCE,
                    ACTIONABLECHARACTERISTIC.SOURCEEVENT,
                    ACTIONABLECHARACTERISTIC.SOURCEURLS,
                    ACTIONABLECHARACTERISTIC.TREATMENT,
                    ACTIONABLECHARACTERISTIC.DRUGCLASSES,
                    ACTIONABLECHARACTERISTIC.APPLICABLECANCERTYPE,
                    ACTIONABLECHARACTERISTIC.APPLICABLEDOID,
                    ACTIONABLECHARACTERISTIC.BLACKLISTCANCERTYPES,
                    ACTIONABLECHARACTERISTIC.LEVEL,
                    ACTIONABLECHARACTERISTIC.DIRECTION,
                    ACTIONABLECHARACTERISTIC.EVIDENCEURLS);
            batch.forEach(entry -> addRecordCharacteristics(timestamp, inserter, entry));
            inserter.execute();
        }

        List<ActionableHLA> actionableHLAs = actionableEvents.hla();
        for (List<ActionableHLA> batch : Iterables.partition(actionableHLAs, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep13 inserter = context.insertInto(ACTIONABLEHLA,
                    ACTIONABLEHLA.MODIFIED,
                    ACTIONABLEHLA.HLATYPE,
                    ACTIONABLEHLA.SOURCE,
                    ACTIONABLEHLA.SOURCEEVENT,
                    ACTIONABLEHLA.SOURCEURLS,
                    ACTIONABLEHLA.TREATMENT,
                    ACTIONABLEHLA.DRUGCLASSES,
                    ACTIONABLEHLA.APPLICABLECANCERTYPE,
                    ACTIONABLEHLA.APPLICABLEDOID,
                    ACTIONABLEHLA.BLACKLISTCANCERTYPES,
                    ACTIONABLEHLA.LEVEL,
                    ACTIONABLEHLA.DIRECTION,
                    ACTIONABLEHLA.EVIDENCEURLS);
            batch.forEach(entry -> addRecordHlas(timestamp, inserter, entry));
            inserter.execute();
        }

        Set<KnownHotspot> knownHotspots = knownEvents.knownHotspots();
        for (List<KnownHotspot> batch : Iterables.partition(knownHotspots, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep9 inserter = context.insertInto(KNOWNHOTSPOT,
                    KNOWNHOTSPOT.MODIFIED,
                    KNOWNHOTSPOT.CHROMOSOME,
                    KNOWNHOTSPOT.POSITION,
                    KNOWNHOTSPOT.REF,
                    KNOWNHOTSPOT.ALT,
                    KNOWNHOTSPOT.INPUTGENE,
                    KNOWNHOTSPOT.INPUTTRANSCRIPT,
                    KNOWNHOTSPOT.INPUTPROTEINANNOTATION,
                    KNOWNHOTSPOT.INPUTSOURCE);
            batch.forEach(entry -> addRecordKnownHotspots(timestamp, inserter, entry));
            inserter.execute();
        }

        Set<KnownCodon> knownCodons = knownEvents.knownCodons();
        for (List<KnownCodon> batch : Iterables.partition(knownCodons, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep9 inserter = context.insertInto(KNOWNCODON,
                    KNOWNCODON.MODIFIED,
                    KNOWNCODON.GENE,
                    KNOWNCODON.TRANSCRIPT,
                    KNOWNCODON.CHROMOSOME,
                    KNOWNCODON.START,
                    KNOWNCODON.END,
                    KNOWNCODON.MUTATIONTYPE,
                    KNOWNCODON.CODONRANK,
                    KNOWNCODON.SOURCES);
            batch.forEach(entry -> addRecordKnownCodons(timestamp, inserter, entry));
            inserter.execute();
        }

        Set<KnownExon> knownExons = knownEvents.knownExons();
        for (List<KnownExon> batch : Iterables.partition(knownExons, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep9 inserter = context.insertInto(KNOWNEXON,
                    KNOWNEXON.MODIFIED,
                    KNOWNEXON.GENE,
                    KNOWNEXON.TRANSCRIPT,
                    KNOWNEXON.CHROMOSOME,
                    KNOWNEXON.START,
                    KNOWNEXON.END,
                    KNOWNEXON.MUTATIONTYPE,
                    KNOWNEXON.EXONRANK,
                    KNOWNEXON.SOURCES);
            batch.forEach(entry -> addRecordKnownExons(timestamp, inserter, entry));
            inserter.execute();
        }

        Set<KnownFusionPair> knownFusionPairs = knownEvents.knownFusionPairs();
        for (List<KnownFusionPair> batch : Iterables.partition(knownFusionPairs, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep8 inserter = context.insertInto(KNOWNFUSIONPAIR,
                    KNOWNFUSIONPAIR.MODIFIED,
                    KNOWNFUSIONPAIR.GENEUP,
                    KNOWNFUSIONPAIR.MINEXONUP,
                    KNOWNFUSIONPAIR.MAXEXONUP,
                    KNOWNFUSIONPAIR.GENEDOWN,
                    KNOWNFUSIONPAIR.MINEXONDOWN,
                    KNOWNFUSIONPAIR.MAXEXONDOWN,
                    KNOWNFUSIONPAIR.SOURCES);
            batch.forEach(entry -> addRecordKnownFusionPairs(timestamp, inserter, entry));
            inserter.execute();
        }

        Set<KnownCopyNumber> knownCopyNumbers = knownEvents.knownCopyNumbers();
        for (List<KnownCopyNumber> batch : Iterables.partition(knownCopyNumbers, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep4 inserter = context.insertInto(KNOWNCOPYNUMBER,
                    KNOWNCOPYNUMBER.MODIFIED,
                    KNOWNCOPYNUMBER.GENE,
                    KNOWNCOPYNUMBER.TYPE,
                    KNOWNCOPYNUMBER.SOURCES);
            batch.forEach(entry -> addRecordKnownCopyNumbers(timestamp, inserter, entry));
            inserter.execute();
        }
    }

    private static void addRecordHotspots(@NotNull Timestamp timestamp, @NotNull InsertValuesStep16 inserter,
            @NotNull ActionableHotspot actionableHotspot) {
        inserter.values(timestamp,
                actionableHotspot.chromosome(),
                actionableHotspot.position(),
                actionableHotspot.ref(),
                actionableHotspot.alt(),
                actionableHotspot.source(),
                actionableHotspot.sourceEvent(),
                ActionableFileFunctions.urlsToString(actionableHotspot.sourceUrls()).isEmpty()
                        ? null
                        : ActionableFileFunctions.urlsToString(actionableHotspot.sourceUrls()),
                actionableHotspot.treatment().treament(),
                ActionableFileFunctions.drugClassesToString(actionableHotspot.treatment().sourceRelevantTreatmentApproaches()).isEmpty()
                        ? null
                        : ActionableFileFunctions.drugClassesToString(actionableHotspot.treatment().sourceRelevantTreatmentApproaches()),
                actionableHotspot.applicableCancerType().name(),
                actionableHotspot.applicableCancerType().doid(),
                CancerTypeFactory.toString(actionableHotspot.blacklistCancerTypes()).isEmpty()
                        ? null
                        : CancerTypeFactory.toString(actionableHotspot.blacklistCancerTypes()),
                actionableHotspot.level(),
                actionableHotspot.direction(),
                ActionableFileFunctions.urlsToString(actionableHotspot.evidenceUrls()).isEmpty()
                        ? null
                        : ActionableFileFunctions.urlsToString(actionableHotspot.evidenceUrls()));
    }

    private static void addRecordRanges(@NotNull Timestamp timestamp, @NotNull InsertValuesStep20 inserter,
            @NotNull ActionableRange actionableRange) {
        inserter.values(timestamp,
                actionableRange.gene(),
                actionableRange.transcript(),
                actionableRange.chromosome(),
                actionableRange.start(),
                actionableRange.end(),
                actionableRange.mutationType(),
                actionableRange.rangeType(),
                actionableRange.rank(),
                actionableRange.source(),
                actionableRange.sourceEvent(),
                ActionableFileFunctions.urlsToString(actionableRange.sourceUrls()).isEmpty()
                        ? null
                        : ActionableFileFunctions.urlsToString(actionableRange.sourceUrls()),
                actionableRange.treatment().treament(),
                ActionableFileFunctions.drugClassesToString(actionableRange.treatment().sourceRelevantTreatmentApproaches()).isEmpty()
                        ? null
                        : ActionableFileFunctions.drugClassesToString(actionableRange.treatment().sourceRelevantTreatmentApproaches()),
                actionableRange.applicableCancerType().name(),
                actionableRange.applicableCancerType().doid(),
                CancerTypeFactory.toString(actionableRange.blacklistCancerTypes()).isEmpty()
                        ? null
                        : CancerTypeFactory.toString(actionableRange.blacklistCancerTypes()),
                actionableRange.level(),
                actionableRange.direction(),
                ActionableFileFunctions.urlsToString(actionableRange.evidenceUrls()).isEmpty()
                        ? null
                        : ActionableFileFunctions.urlsToString(actionableRange.evidenceUrls()));
    }

    private static void addRecordGenes(@NotNull Timestamp timestamp, @NotNull InsertValuesStep14 inserter,
            @NotNull ActionableGene actionableGene) {
        inserter.values(timestamp,
                actionableGene.gene(),
                actionableGene.event(),
                actionableGene.source(),
                actionableGene.sourceEvent(),
                ActionableFileFunctions.urlsToString(actionableGene.sourceUrls()).isEmpty()
                        ? null
                        : ActionableFileFunctions.urlsToString(actionableGene.sourceUrls()),
                actionableGene.treatment().treament(),
                ActionableFileFunctions.drugClassesToString(actionableGene.treatment().sourceRelevantTreatmentApproaches()).isEmpty()
                        ? null
                        : ActionableFileFunctions.drugClassesToString(actionableGene.treatment().sourceRelevantTreatmentApproaches()),
                actionableGene.applicableCancerType().name(),
                actionableGene.applicableCancerType().doid(),
                CancerTypeFactory.toString(actionableGene.blacklistCancerTypes()).isEmpty()
                        ? null
                        : CancerTypeFactory.toString(actionableGene.blacklistCancerTypes()),
                actionableGene.level(),
                actionableGene.direction(),
                ActionableFileFunctions.urlsToString(actionableGene.evidenceUrls()).isEmpty()
                        ? null
                        : ActionableFileFunctions.urlsToString(actionableGene.evidenceUrls()));
    }

    private static void addRecordFusions(@NotNull Timestamp timestamp, @NotNull InsertValuesStep18 inserter,
            @NotNull ActionableFusion actionableFusion) {
        inserter.values(timestamp,
                actionableFusion.geneUp(),
                actionableFusion.minExonUp(),
                actionableFusion.maxExonUp(),
                actionableFusion.geneDown(),
                actionableFusion.minExonDown(),
                actionableFusion.maxExonDown(),
                actionableFusion.source(),
                actionableFusion.sourceEvent(),
                ActionableFileFunctions.urlsToString(actionableFusion.sourceUrls()).isEmpty()
                        ? null
                        : ActionableFileFunctions.urlsToString(actionableFusion.sourceUrls()),
                actionableFusion.treatment().treament(),
                ActionableFileFunctions.drugClassesToString(actionableFusion.treatment().sourceRelevantTreatmentApproaches()).isEmpty()
                        ? null
                        : ActionableFileFunctions.drugClassesToString(actionableFusion.treatment().sourceRelevantTreatmentApproaches()),
                actionableFusion.applicableCancerType().name(),
                actionableFusion.applicableCancerType().doid(),
                CancerTypeFactory.toString(actionableFusion.blacklistCancerTypes()).isEmpty()
                        ? null
                        : CancerTypeFactory.toString(actionableFusion.blacklistCancerTypes()),
                actionableFusion.level(),
                actionableFusion.direction(),
                ActionableFileFunctions.urlsToString(actionableFusion.evidenceUrls()).isEmpty()
                        ? null
                        : ActionableFileFunctions.urlsToString(actionableFusion.evidenceUrls()));
    }

    private static void addRecordCharacteristics(@NotNull Timestamp timestamp, @NotNull InsertValuesStep15 inserter,
            @NotNull ActionableCharacteristic actionableCharacteristic) {
        inserter.values(timestamp,
                actionableCharacteristic.name(),
                actionableCharacteristic.comparator(),
                actionableCharacteristic.cutoff(),
                actionableCharacteristic.source(),
                actionableCharacteristic.sourceEvent(),
                ActionableFileFunctions.urlsToString(actionableCharacteristic.sourceUrls()).isEmpty()
                        ? null
                        : ActionableFileFunctions.urlsToString(actionableCharacteristic.sourceUrls()),
                actionableCharacteristic.treatment().treament(),
                ActionableFileFunctions.drugClassesToString(actionableCharacteristic.treatment().sourceRelevantTreatmentApproaches())
                        .isEmpty()
                        ? null
                        : ActionableFileFunctions.drugClassesToString(actionableCharacteristic.treatment()
                                .sourceRelevantTreatmentApproaches()),
                actionableCharacteristic.applicableCancerType().name(),
                actionableCharacteristic.applicableCancerType().doid(),
                CancerTypeFactory.toString(actionableCharacteristic.blacklistCancerTypes()).isEmpty()
                        ? null
                        : CancerTypeFactory.toString(actionableCharacteristic.blacklistCancerTypes()),
                actionableCharacteristic.level(),
                actionableCharacteristic.direction(),
                ActionableFileFunctions.urlsToString(actionableCharacteristic.evidenceUrls()).isEmpty()
                        ? null
                        : ActionableFileFunctions.urlsToString(actionableCharacteristic.evidenceUrls()));
    }

    private static void addRecordHlas(@NotNull Timestamp timestamp, @NotNull InsertValuesStep13 inserter,
            @NotNull ActionableHLA actionableHLA) {
        inserter.values(timestamp,
                actionableHLA.hlaType(),
                actionableHLA.source(),
                actionableHLA.sourceEvent(),
                ActionableFileFunctions.urlsToString(actionableHLA.sourceUrls()).isEmpty()
                        ? null
                        : ActionableFileFunctions.urlsToString(actionableHLA.sourceUrls()),
                actionableHLA.treatment().treament(),
                ActionableFileFunctions.drugClassesToString(actionableHLA.treatment().sourceRelevantTreatmentApproaches()).isEmpty()
                        ? null
                        : ActionableFileFunctions.drugClassesToString(actionableHLA.treatment().sourceRelevantTreatmentApproaches()),
                actionableHLA.applicableCancerType().name(),
                actionableHLA.applicableCancerType().doid(),
                CancerTypeFactory.toString(actionableHLA.blacklistCancerTypes()).isEmpty()
                        ? null
                        : CancerTypeFactory.toString(actionableHLA.blacklistCancerTypes()),
                actionableHLA.level(),
                actionableHLA.direction(),
                ActionableFileFunctions.urlsToString(actionableHLA.evidenceUrls()).isEmpty()
                        ? null
                        : ActionableFileFunctions.urlsToString(actionableHLA.evidenceUrls()));
    }

    private static void addRecordKnownHotspots(@NotNull Timestamp timestamp, @NotNull InsertValuesStep9 inserter,
            @NotNull KnownHotspot knownHotspot) {
        inserter.values(timestamp,
                knownHotspot.chromosome(),
                knownHotspot.position(),
                knownHotspot.ref(),
                knownHotspot.alt(),
                knownHotspot.gene(),
                knownHotspot.transcript(),
                knownHotspot.proteinAnnotation(),
                Knowledgebase.toCommaSeparatedSourceString(knownHotspot.sources()));
    }

    private static void addRecordKnownCodons(@NotNull Timestamp timestamp, @NotNull InsertValuesStep9 inserter,
            @NotNull KnownCodon knownCodon) {
        inserter.values(timestamp,
                knownCodon.annotation().gene(),
                knownCodon.annotation().transcript(),
                knownCodon.annotation().chromosome(),
                knownCodon.annotation().start(),
                knownCodon.annotation().end(),
                knownCodon.annotation().mutationType(),
                knownCodon.annotation().rank(),
                Knowledgebase.toCommaSeparatedSourceString(knownCodon.sources()));
    }

    private static void addRecordKnownExons(@NotNull Timestamp timestamp, @NotNull InsertValuesStep9 inserter,
            @NotNull KnownExon knownExon) {
        inserter.values(timestamp,
                knownExon.annotation().gene(),
                knownExon.annotation().transcript(),
                knownExon.annotation().chromosome(),
                knownExon.annotation().start(),
                knownExon.annotation().end(),
                knownExon.annotation().mutationType(),
                knownExon.annotation().rank(),
                Knowledgebase.toCommaSeparatedSourceString(knownExon.sources()));
    }

    private static void addRecordKnownFusionPairs(@NotNull Timestamp timestamp, @NotNull InsertValuesStep8 inserter,
            @NotNull KnownFusionPair knownFusionPairs) {
        inserter.values(timestamp,
                knownFusionPairs.geneUp(),
                knownFusionPairs.minExonUp(),
                knownFusionPairs.maxExonUp(),
                knownFusionPairs.geneDown(),
                knownFusionPairs.minExonDown(),
                knownFusionPairs.maxExonDown(),
                Knowledgebase.toCommaSeparatedSourceString(knownFusionPairs.sources()));
    }

    private static void addRecordKnownCopyNumbers(@NotNull Timestamp timestamp, @NotNull InsertValuesStep4 inserter,
            @NotNull KnownCopyNumber knownCopyNumber) {
        inserter.values(timestamp,
                knownCopyNumber.gene(),
                knownCopyNumber.type(),
                Knowledgebase.toCommaSeparatedSourceString(knownCopyNumber.sources()));
    }
}