package com.hartwig.hmftools.patientdb;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Set;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.driver.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.driver.LikelihoodMethod;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.FittedPurityScore;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineDetectionMethod;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurityScore;
import com.hartwig.hmftools.common.purple.ImmutablePurityContext;
import com.hartwig.hmftools.common.purple.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.common.purple.MicrosatelliteStatus;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.purple.RunMode;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.purple.TumorMutationalStatus;
import com.hartwig.hmftools.common.sv.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.CopynumberRecord;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.DrivercatalogRecord;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.GenecopynumberRecord;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.GermlinedeletionRecord;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.PurityRecord;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.StructuralvariantRecord;

import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class LoadPurpleDataTest extends DatabaseTestBase
{
    @Test
    public void canWritePurityContext()
    {
        PurpleQC qc = ImmutablePurpleQC.builder()
                .status(List.of(PurpleQCStatus.PASS))
                .method(FittedPurityMethod.SOMATIC)
                .copyNumberSegments(0)
                .unsupportedCopyNumberSegments(0)
                .deletedGenes(0)
                .purity(0)
                .contamination(0)
                .amberMeanDepth(0)
                .cobaltGender(Gender.FEMALE)
                .amberGender(Gender.FEMALE)
                .germlineAberrations(Set.of(GermlineAberration.NONE))
                .lohPercent(0)
                .chimerismPercentage(0)
                .tincLevel(0)
                .build();

        FittedPurity fittedPurity = ImmutableFittedPurity.builder()
                .purity(0)
                .normFactor(0)
                .score(0)
                .diploidProportion(0)
                .ploidy(0)
                .somaticPenalty(0)
                .build();

        FittedPurityScore score = ImmutableFittedPurityScore.builder()
                .minPurity(0)
                .maxPurity(0)
                .minPloidy(0)
                .maxPloidy(0)
                .minDiploidProportion(0)
                .maxDiploidProportion(0)
                .build();

        PurityContext purityContext = ImmutablePurityContext.builder()
                .qc(qc)
                .bestFit(fittedPurity)
                .score(score)
                .gender(Gender.FEMALE)
                .method(FittedPurityMethod.SOMATIC)
                .runMode(RunMode.TUMOR_GERMLINE)
                .targeted(false)
                .polyClonalProportion(0)
                .wholeGenomeDuplication(false)
                .microsatelliteIndelsPerMb(0)
                .microsatelliteStatus(MicrosatelliteStatus.MSS)
                .tumorMutationalLoad(0)
                .tumorMutationalLoadStatus(TumorMutationalStatus.LOW)
                .tumorMutationalBurdenPerMb(0)
                .tumorMutationalBurdenStatus(TumorMutationalStatus.LOW)
                .svTumorMutationalBurden(0)
                .build();

        databaseAccess.writePurity(TEST_SAMPLE_ID, purityContext, purityContext.qc());

        List<PurityRecord> purityRecords = fetchTable(Tables.PURITY, PurityRecord.class);
        assertEquals(1, purityRecords.size());
    }

    @Test
    public void canWriteCopyNumbers()
    {
        PurpleCopyNumber copyNumber = ImmutablePurpleCopyNumber.builder()
                .chromosome("chr1")
                .start(1)
                .end(1)
                .bafCount(0)
                .averageActualBAF(0)
                .averageObservedBAF(0)
                .averageTumorCopyNumber(0)
                .depthWindowCount(0)
                .segmentStartSupport(SegmentSupport.NONE)
                .segmentEndSupport(SegmentSupport.NONE)
                .method(CopyNumberMethod.UNKNOWN)
                .gcContent(0)
                .minStart(0)
                .maxStart(0)
                .build();

        databaseAccess.writeCopynumbers(TEST_SAMPLE_ID, List.of(copyNumber));

        List<CopynumberRecord> copyNumberRecords = fetchTable(Tables.COPYNUMBER, CopynumberRecord.class);
        assertEquals(1, copyNumberRecords.size());
    }

    @Test
    public void canWriteGeneCopyNumber()
    {
        GeneCopyNumber geneCopyNumber = new GeneCopyNumber(
                "chr1",
                1,
                1,
                "GENE1",
                "ENST0",
                true,
                "1q1.1",
                1,
                0,
                0,
                0,
                0,
                1,
                1,
                0,
                0,
                SegmentSupport.NONE,
                SegmentSupport.NONE,
                CopyNumberMethod.UNKNOWN,
                0
        );

        databaseAccess.writeGeneCopyNumbers(TEST_SAMPLE_ID, List.of(geneCopyNumber));

        List<GenecopynumberRecord> geneCopyNumberRecords = fetchTable(Tables.GENECOPYNUMBER, GenecopynumberRecord.class);
        assertEquals(1, geneCopyNumberRecords.size());
    }

    @Test
    public void canWriteDriverCatalog()
    {
        List<DriverCatalog> drivers = Stream.of(DriverType.values()).map(x -> createDriverWithType(x)).toList();

        databaseAccess.writePurpleDriverCatalog(TEST_SAMPLE_ID, drivers, null);

        List<DrivercatalogRecord> driverCatalogRecords = fetchTable(Tables.DRIVERCATALOG, DrivercatalogRecord.class);
        assertEquals(DriverType.DRIVERS_PURPLE_SOMATIC.size(), driverCatalogRecords.size());
    }

    private static DriverCatalog createDriverWithType(DriverType driverType)
    {
        return ImmutableDriverCatalog.builder()
                .chromosome("chr1")
                .chromosomeBand("1q1.1")
                .gene("GENE1")
                .transcript("ENST0")
                .isCanonical(true)
                .category(DriverCategory.TSG)
                .driver(driverType)
                .likelihoodMethod(LikelihoodMethod.NONE)
                .reportedStatus(ReportedStatus.NONE)
                .driverLikelihood(0)
                .missense(0)
                .nonsense(0)
                .splice(0)
                .inframe(0)
                .frameshift(0)
                .biallelic(false)
                .minCopyNumber(0)
                .maxCopyNumber(0)
                .build();
    }

    @Test
    public void canWriteGermlineDeletions()
    {
        GermlineDeletion germlineDeletion = new GermlineDeletion(
                "GENE1",
                "chr1",
                "1q1.1",
                1,
                1,
                1,
                1,
                1,
                GermlineDetectionMethod.SEGMENT,
                GermlineStatus.UNKNOWN,
                GermlineStatus.UNKNOWN,
                0,
                0,
                "PASS",
                0,
                ReportedStatus.REPORTED
        );

        databaseAccess.writeGermlineDeletions(TEST_SAMPLE_ID, List.of(germlineDeletion));

        List<GermlinedeletionRecord> germlineDeletionRecords = fetchTable(Tables.GERMLINEDELETION, GermlinedeletionRecord.class);
        assertEquals(1, germlineDeletionRecords.size());
    }

    @Test
    public void canWriteSmallVariants()
    {
        // TODO
    }

    @Test
    public void canWriteStructuralVariants()
    {
        ImmutableStructuralVariantData structuralVariant = ImmutableStructuralVariantData.builder()
                .id(0)
                .vcfIdStart("0")
                .vcfIdEnd("1")
                .startChromosome("chr1")
                .endChromosome("chr1")
                .startPosition(1)
                .endPosition(2)
                .startOrientation((byte) 1)
                .endOrientation((byte) 1)
                .startHomologySequence("A")
                .endHomologySequence("A")
                .junctionCopyNumber(0)
                .startAF(0)
                .endAF(0)
                .adjustedStartAF(0)
                .adjustedEndAF(0)
                .adjustedStartCopyNumber(0)
                .adjustedEndCopyNumber(0)
                .adjustedStartCopyNumberChange(0)
                .adjustedEndCopyNumberChange(0)
                .insertSequence("A")
                .type(StructuralVariantType.BND)
                .filter("PASS")
                .qualityScore(0)
                .event("")
                .startTumorVariantFragmentCount(0)
                .startTumorReferenceFragmentCount(0)
                .startNormalVariantFragmentCount(0)
                .startNormalReferenceFragmentCount(0)
                .endTumorVariantFragmentCount(0)
                .endTumorReferenceFragmentCount(0)
                .endNormalVariantFragmentCount(0)
                .endNormalReferenceFragmentCount(0)
                .startIntervalOffsetStart(0)
                .startIntervalOffsetEnd(0)
                .endIntervalOffsetStart(0)
                .endIntervalOffsetEnd(0)
                .inexactHomologyOffsetStart(0)
                .inexactHomologyOffsetEnd(0)
                .startLinkedBy("0")
                .endLinkedBy("1")
                .insertSequenceAlignments("")
                .insertSequenceRepeatClass("")
                .insertSequenceRepeatType("")
                .insertSequenceRepeatOrientation((byte) 0 )
                .insertSequenceRepeatCoverage(0)
                .startAnchoringSupportDistance(0)
                .endAnchoringSupportDistance(0)
                .ponCount(0)
                .build();

        databaseAccess.writeStructuralVariants(TEST_SAMPLE_ID, List.of(structuralVariant));

        List<StructuralvariantRecord> structuralVariantRecords = fetchTable(Tables.STRUCTURALVARIANT, StructuralvariantRecord.class);
        assertEquals(1, structuralVariantRecords.size());
    }
}
