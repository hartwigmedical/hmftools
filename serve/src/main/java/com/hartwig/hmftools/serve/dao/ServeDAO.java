package com.hartwig.hmftools.serve.dao;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.serve.actionability.ActionableEvents;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.immuno.ActionableHLA;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.database.tables.Knownfusionpairs;
import com.hartwig.hmftools.serve.extraction.KnownEvents;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodon;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPair;

import static com.hartwig.hmftools.serve.database.tables.Actionablehotspots.ACTIONABLEHOTSPOTS;
import static com.hartwig.hmftools.serve.database.tables.Actionableranges.ACTIONABLERANGES;
import static com.hartwig.hmftools.serve.database.tables.Actionablegenes.ACTIONABLEGENES;
import static com.hartwig.hmftools.serve.database.tables.Actionablefusions.ACTIONABLEFUSIONS;
import static com.hartwig.hmftools.serve.database.tables.Actionablecharacteristics.ACTIONABLECHARACTERISTICS;
import static com.hartwig.hmftools.serve.database.tables.Actionablehla.ACTIONABLEHLA;
import static com.hartwig.hmftools.serve.database.tables.Knowncodons.KNOWNCODONS;
import static com.hartwig.hmftools.serve.database.tables.Knownexons.KNOWNEXONS;
import static com.hartwig.hmftools.serve.database.tables.Knowncopynumbers.KNOWNCOPYNUMBERS;
import static com.hartwig.hmftools.serve.database.tables.Knownfusionpairs.KNOWNFUSIONPAIRS;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep13;
import org.jooq.InsertValuesStep14;
import org.jooq.InsertValuesStep15;
import org.jooq.InsertValuesStep16;
import org.jooq.InsertValuesStep18;
import org.jooq.InsertValuesStep20;
import org.jooq.InsertValuesStep4;
import org.jooq.InsertValuesStep7;
import org.jooq.InsertValuesStep9;

public class ServeDAO {

    @NotNull
    private final DSLContext context;

    ServeDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void write(@NotNull ActionableEvents actionableEvents, KnownEvents knownEvents) {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        List<ActionableHotspot> actionableHotspots = actionableEvents.hotspots();
        for (List<ActionableHotspot> batch : Iterables.partition(actionableHotspots, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep16 inserter = context.insertInto(ACTIONABLEHOTSPOTS,
                    ACTIONABLEHOTSPOTS.MODIFIED,
                    ACTIONABLEHOTSPOTS.CHROMOSOME,
                    ACTIONABLEHOTSPOTS.POSITION,
                    ACTIONABLEHOTSPOTS.REF,
                    ACTIONABLEHOTSPOTS.ALT,
                    ACTIONABLEHOTSPOTS.SOURCE,
                    ACTIONABLEHOTSPOTS.SOURCEEVENT,
                    ACTIONABLEHOTSPOTS.SOURCEURLS,
                    ACTIONABLEHOTSPOTS.TREATMENT,
                    ACTIONABLEHOTSPOTS.DRUGCLASSES,
                    ACTIONABLEHOTSPOTS.APPLICABLECANCERTYPE,
                    ACTIONABLEHOTSPOTS.APPLICABLEDOID,
                    ACTIONABLEHOTSPOTS.BLACKLISTCANCERTYPES,
                    ACTIONABLEHOTSPOTS.LEVEL,
                    ACTIONABLEHOTSPOTS.DIRECTION,
                    ACTIONABLEHOTSPOTS.EVIDENCEURLS);
            batch.forEach(entry -> addRecordHotspots(timestamp, inserter, entry));
            inserter.execute();
        }

        List<ActionableRange> actionableRanges = actionableEvents.ranges();
        for (List<ActionableRange> batch : Iterables.partition(actionableRanges, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep20 inserter = context.insertInto(ACTIONABLERANGES,
                    ACTIONABLERANGES.MODIFIED,
                    ACTIONABLERANGES.GENE,
                    ACTIONABLERANGES.TRANSCRIPT,
                    ACTIONABLERANGES.CHROMOSOME,
                    ACTIONABLERANGES.START,
                    ACTIONABLERANGES.END,
                    ACTIONABLERANGES.MUTATIONTYPE,
                    ACTIONABLERANGES.RANGETYPE,
                    ACTIONABLERANGES.RANK,
                    ACTIONABLERANGES.SOURCE,
                    ACTIONABLERANGES.SOURCEEVENT,
                    ACTIONABLERANGES.SOURCEURLS,
                    ACTIONABLERANGES.TREATMENT,
                    ACTIONABLERANGES.DRUGCLASSES,
                    ACTIONABLERANGES.APPLICABLECANCERTYPE,
                    ACTIONABLERANGES.APPLICABLEDOID,
                    ACTIONABLERANGES.BLACKLISTCANCERTYPES,
                    ACTIONABLERANGES.LEVEL,
                    ACTIONABLERANGES.DIRECTION,
                    ACTIONABLERANGES.EVIDENCEURLS);
            batch.forEach(entry -> addRecordRanges(timestamp, inserter, entry));
            inserter.execute();
        }

        List<ActionableGene> actionableGenes = actionableEvents.genes();
        for (List<ActionableGene> batch : Iterables.partition(actionableGenes, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep14 inserter = context.insertInto(ACTIONABLEGENES,
                    ACTIONABLEGENES.MODIFIED,
                    ACTIONABLEGENES.GENE,
                    ACTIONABLEGENES.EVENT,
                    ACTIONABLEGENES.SOURCE,
                    ACTIONABLEGENES.SOURCEEVENT,
                    ACTIONABLEGENES.SOURCEURLS,
                    ACTIONABLEGENES.TREATMENT,
                    ACTIONABLEGENES.DRUGCLASSES,
                    ACTIONABLEGENES.APPLICABLECANCERTYPE,
                    ACTIONABLEGENES.APPLICABLEDOID,
                    ACTIONABLEGENES.BLACKLISTCANCERTYPES,
                    ACTIONABLEGENES.LEVEL,
                    ACTIONABLEGENES.DIRECTION,
                    ACTIONABLEGENES.EVIDENCEURLS);
            batch.forEach(entry -> addRecordGenes(timestamp, inserter, entry));
            inserter.execute();
        }

        List<ActionableFusion> actionableFusions = actionableEvents.fusions();
        for (List<ActionableFusion> batch : Iterables.partition(actionableFusions, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep18 inserter = context.insertInto(ACTIONABLEFUSIONS,
                    ACTIONABLEFUSIONS.MODIFIED,
                    ACTIONABLEFUSIONS.GENEUP,
                    ACTIONABLEFUSIONS.MINEXONUP,
                    ACTIONABLEFUSIONS.MAXEXONUP,
                    ACTIONABLEFUSIONS.GENEDOWN,
                    ACTIONABLEFUSIONS.MINEXONDOWN,
                    ACTIONABLEFUSIONS.MAXEXONDOWN,
                    ACTIONABLEFUSIONS.SOURCE,
                    ACTIONABLEFUSIONS.SOURCEEVENT,
                    ACTIONABLEFUSIONS.SOURCEURLS,
                    ACTIONABLEFUSIONS.TREATMENT,
                    ACTIONABLEFUSIONS.DRUGCLASSES,
                    ACTIONABLEFUSIONS.APPLICABLECANCERTYPE,
                    ACTIONABLEFUSIONS.APPLICABLEDOID,
                    ACTIONABLEFUSIONS.BLACKLISTCANCERTYPES,
                    ACTIONABLEFUSIONS.LEVEL,
                    ACTIONABLEFUSIONS.DIRECTION,
                    ACTIONABLEFUSIONS.EVIDENCEURLS);
            batch.forEach(entry -> addRecordFusions(timestamp, inserter, entry));
            inserter.execute();
        }

        List<ActionableCharacteristic> actionableCharacteristics = actionableEvents.characteristics();
        for (List<ActionableCharacteristic> batch : Iterables.partition(actionableCharacteristics, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep15 inserter = context.insertInto(ACTIONABLECHARACTERISTICS,
                    ACTIONABLECHARACTERISTICS.MODIFIED,
                    ACTIONABLECHARACTERISTICS.NAME,
                    ACTIONABLECHARACTERISTICS.COMPARATOR,
                    ACTIONABLECHARACTERISTICS.CUTOFF,
                    ACTIONABLECHARACTERISTICS.SOURCE,
                    ACTIONABLECHARACTERISTICS.SOURCEEVENT,
                    ACTIONABLECHARACTERISTICS.SOURCEURLS,
                    ACTIONABLECHARACTERISTICS.TREATMENT,
                    ACTIONABLECHARACTERISTICS.DRUGCLASSES,
                    ACTIONABLECHARACTERISTICS.APPLICABLECANCERTYPE,
                    ACTIONABLECHARACTERISTICS.APPLICABLEDOID,
                    ACTIONABLECHARACTERISTICS.BLACKLISTCANCERTYPES,
                    ACTIONABLECHARACTERISTICS.LEVEL,
                    ACTIONABLECHARACTERISTICS.DIRECTION,
                    ACTIONABLECHARACTERISTICS.EVIDENCEURLS);
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

        Set<KnownCodon> knownCodons = knownEvents.knownCodons();
        for (List<KnownCodon> batch : Iterables.partition(knownCodons, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep9 inserter = context.insertInto(KNOWNCODONS,
                    KNOWNCODONS.MODIFIED,
                    KNOWNCODONS.GENE,
                    KNOWNCODONS.TRANSCRIPT,
                    KNOWNCODONS.CHROMOSOME,
                    KNOWNCODONS.START,
                    KNOWNCODONS.END,
                    KNOWNCODONS.MUTATIONTYPE,
                    KNOWNCODONS.CODONRANK,
                    KNOWNCODONS.SOURCES);
            batch.forEach(entry -> addRecordKnownCodons(timestamp, inserter, entry));
            inserter.execute();
        }

        Set<KnownExon> knownExons = knownEvents.knownExons();
        for (List<KnownExon> batch : Iterables.partition(knownExons, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep9 inserter = context.insertInto(KNOWNEXONS,
                    KNOWNEXONS.MODIFIED,
                    KNOWNEXONS.GENE,
                    KNOWNEXONS.TRANSCRIPT,
                    KNOWNEXONS.CHROMOSOME,
                    KNOWNEXONS.START,
                    KNOWNEXONS.END,
                    KNOWNEXONS.MUTATIONTYPE,
                    KNOWNEXONS.EXONRANK,
                    KNOWNEXONS.SOURCES);
            batch.forEach(entry -> addRecordKnownExons(timestamp, inserter, entry));
            inserter.execute();
        }

        Set<KnownFusionPair> knownFusionPairs = knownEvents.knownFusionPairs();
        for (List<KnownFusionPair> batch : Iterables.partition(knownFusionPairs, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep7 inserter = context.insertInto(KNOWNFUSIONPAIRS,
                    KNOWNFUSIONPAIRS.MODIFIED,
                    KNOWNFUSIONPAIRS.GENEUP,
                    KNOWNFUSIONPAIRS.MINEXONUP,
                    KNOWNFUSIONPAIRS.MAXEXONUP,
                    KNOWNFUSIONPAIRS.GENEDOWN,
                    KNOWNFUSIONPAIRS.MINEXONDOWN,
                    KNOWNFUSIONPAIRS.MAXEXONDOWN);
            batch.forEach(entry -> addRecordKnownFusionPairs(timestamp, inserter, entry));
            inserter.execute();
        }

        Set<KnownCopyNumber> knownCopyNumbers = knownEvents.knownCopyNumbers();
        for (List<KnownCopyNumber> batch : Iterables.partition(knownCopyNumbers, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep4 inserter = context.insertInto(KNOWNCOPYNUMBERS,
                    KNOWNCOPYNUMBERS.MODIFIED,
                    KNOWNCOPYNUMBERS.GENE,
                    KNOWNCOPYNUMBERS.TYPE,
                    KNOWNCOPYNUMBERS.SOURCES);
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
                actionableHotspot.sourceUrls(),
                actionableHotspot.treatment().treament(),
                actionableHotspot.treatment().drugClasses(),
                actionableHotspot.applicableCancerType().name(),
                actionableHotspot.applicableCancerType().doid(),
                actionableHotspot.blacklistCancerTypes(),
                actionableHotspot.level(),
                actionableHotspot.direction(),
                actionableHotspot.evidenceUrls());
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
                actionableRange.sourceUrls(),
                actionableRange.treatment().treament(),
                actionableRange.treatment().drugClasses(),
                actionableRange.applicableCancerType().name(),
                actionableRange.applicableCancerType().doid(),
                actionableRange.blacklistCancerTypes(),
                actionableRange.level(),
                actionableRange.direction(),
                actionableRange.evidenceUrls());
    }

    private static void addRecordGenes(@NotNull Timestamp timestamp, @NotNull InsertValuesStep14 inserter,
            @NotNull ActionableGene actionableGene) {
        inserter.values(timestamp,
                actionableGene.gene(),
                actionableGene.event(),
                actionableGene.source(),
                actionableGene.sourceEvent(),
                actionableGene.sourceUrls(),
                actionableGene.treatment().treament(),
                actionableGene.treatment().drugClasses(),
                actionableGene.applicableCancerType().name(),
                actionableGene.applicableCancerType().doid(),
                actionableGene.blacklistCancerTypes(),
                actionableGene.level(),
                actionableGene.direction(),
                actionableGene.evidenceUrls());
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
                actionableFusion.sourceUrls(),
                actionableFusion.treatment().treament(),
                actionableFusion.treatment().drugClasses(),
                actionableFusion.applicableCancerType().name(),
                actionableFusion.applicableCancerType().doid(),
                actionableFusion.blacklistCancerTypes(),
                actionableFusion.level(),
                actionableFusion.direction(),
                actionableFusion.evidenceUrls());
    }

    private static void addRecordCharacteristics(@NotNull Timestamp timestamp, @NotNull InsertValuesStep15 inserter,
            @NotNull ActionableCharacteristic actionableCharacteristic) {
        inserter.values(timestamp,
                actionableCharacteristic.name(),
                actionableCharacteristic.comparator(),
                actionableCharacteristic.cutoff(),
                actionableCharacteristic.source(),
                actionableCharacteristic.sourceEvent(),
                actionableCharacteristic.sourceUrls(),
                actionableCharacteristic.treatment().treament(),
                actionableCharacteristic.treatment().drugClasses(),
                actionableCharacteristic.applicableCancerType().name(),
                actionableCharacteristic.applicableCancerType().doid(),
                actionableCharacteristic.blacklistCancerTypes(),
                actionableCharacteristic.level(),
                actionableCharacteristic.direction(),
                actionableCharacteristic.evidenceUrls());
    }

    private static void addRecordHlas(@NotNull Timestamp timestamp, @NotNull InsertValuesStep13 inserter,
            @NotNull ActionableHLA actionableHLA) {
        inserter.values(timestamp,
                actionableHLA.hlaType(),
                actionableHLA.source(),
                actionableHLA.sourceEvent(),
                actionableHLA.sourceUrls(),
                actionableHLA.treatment().treament(),
                actionableHLA.treatment().drugClasses(),
                actionableHLA.applicableCancerType().name(),
                actionableHLA.applicableCancerType().doid(),
                actionableHLA.blacklistCancerTypes(),
                actionableHLA.level(),
                actionableHLA.direction(),
                actionableHLA.evidenceUrls());
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
                knownCodon.sources());
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
                knownExon.sources());
    }

    private static void addRecordKnownFusionPairs(@NotNull Timestamp timestamp, @NotNull InsertValuesStep7 inserter,
            @NotNull KnownFusionPair knownFusionPairs) {
        inserter.values(timestamp,
                knownFusionPairs.geneUp(),
                knownFusionPairs.minExonUp(),
                knownFusionPairs.maxExonUp(),
                knownFusionPairs.geneDown(),
                knownFusionPairs.minExonDown(),
                knownFusionPairs.maxExonDown());
    }

    private static void addRecordKnownCopyNumbers(@NotNull Timestamp timestamp, @NotNull InsertValuesStep4 inserter,
            @NotNull KnownCopyNumber knownCopyNumber) {
        inserter.values(timestamp, knownCopyNumber.gene(), knownCopyNumber.type(), knownCopyNumber.sources());
    }
}