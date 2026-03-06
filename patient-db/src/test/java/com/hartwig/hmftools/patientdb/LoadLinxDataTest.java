package com.hartwig.hmftools.patientdb;

import static org.junit.Assert.assertEquals;

import java.util.EnumSet;
import java.util.List;

import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.driver.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.driver.LikelihoodMethod;
import com.hartwig.hmftools.common.gene.TranscriptCodingType;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;
import com.hartwig.hmftools.common.linx.DriverEventType;
import com.hartwig.hmftools.common.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.linx.FusionPhasedType;
import com.hartwig.hmftools.common.linx.FusionReportableReason;
import com.hartwig.hmftools.common.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.common.linx.ImmutableLinxCluster;
import com.hartwig.hmftools.common.linx.ImmutableLinxDriver;
import com.hartwig.hmftools.common.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.linx.ImmutableLinxLink;
import com.hartwig.hmftools.common.linx.ImmutableLinxSvAnnotation;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxCluster;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.linx.LinxLink;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.DrivercatalogRecord;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.SvannotationRecord;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.SvbreakendRecord;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.SvclusterRecord;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.SvdriverRecord;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.SvfusionRecord;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.SvlinkRecord;

import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class LoadLinxDataTest extends DatabaseTestBase
{
    @Test
    public void canWriteLinxAnnotations()
    {
        LinxSvAnnotation annotation = ImmutableLinxSvAnnotation.builder()
                .svId(0)
                .vcfIdStart("0")
                .vcfIdEnd("1")
                .coordsStart("chr1:1")
                .coordsEnd("chr1:2")
                .type(StructuralVariantType.BND)
                .clusterId(0)
                .clusterReason("")
                .fragileSiteStart(false)
                .fragileSiteEnd(false)
                .isFoldback(false)
                .lineTypeStart("")
                .lineTypeEnd("")
                .junctionCopyNumberMin(0)
                .junctionCopyNumberMax(0)
                .geneStart("")
                .geneEnd("")
                .localTopologyIdStart(0)
                .localTopologyIdEnd(0)
                .localTopologyStart("")
                .localTopologyEnd("")
                .localTICountStart(0)
                .localTICountEnd(0)
                .build();

        databaseAccess.writeSvLinxData(TEST_SAMPLE_ID, List.of(annotation));
        List<SvannotationRecord> records = fetchTable(Tables.SVANNOTATION, SvannotationRecord.class);
        assertEquals(1, records.size());
    }

    @Test
    public void canWriteClusters()
    {
        LinxCluster cluster = ImmutableLinxCluster.builder()
                .clusterId(0)
                .category("COMPLEX")
                .synthetic(false)
                .resolvedType("COMPLEX")
                .clusterCount(3)
                .clusterDesc("DEL=1_INV=2")
                .build();

        databaseAccess.writeSvClusters(TEST_SAMPLE_ID, List.of(cluster));
        List<SvclusterRecord> records = fetchTable(Tables.SVCLUSTER, SvclusterRecord.class);
        assertEquals(1, records.size());
    }

    @Test
    public void canWriteLinks()
    {
        LinxLink link = ImmutableLinxLink.builder()
                .clusterId(0)
                .chainId(0)
                .chainIndex("0")
                .chainCount(1)
                .lowerSvId(0)
                .upperSvId(1)
                .lowerBreakendIsStart(true)
                .upperBreakendIsStart(false)
                .chromosome("chr1")
                .arm("P")
                .assembled(false)
                .traversedSVCount(1)
                .length(1)
                .junctionCopyNumber(1)
                .junctionCopyNumberUncertainty(0.1)
                .pseudogeneInfo("")
                .ecDna(false)
                .build();

        databaseAccess.writeSvLinks(TEST_SAMPLE_ID, List.of(link));
        List<SvlinkRecord> records = fetchTable(Tables.SVLINK, SvlinkRecord.class);
        assertEquals(1, records.size());
    }

    @Test
    public void canWriteBreakendsAndFusions()
    {
        final int FIVE_PRIME_BREAKEND_ID = 0;
        final int THREE_PRIME_BREAKEND_ID = 1;

        LinxBreakend breakend5 = createBreakendWithId(FIVE_PRIME_BREAKEND_ID);
        LinxBreakend breakend3 = createBreakendWithId(THREE_PRIME_BREAKEND_ID);

        LinxFusion fusion = ImmutableLinxFusion.builder()
                .fivePrimeBreakendId(FIVE_PRIME_BREAKEND_ID)
                .threePrimeBreakendId(THREE_PRIME_BREAKEND_ID)
                .name("FUSION")
                .reported(true)
                .reportedType("FUSION")
                .reportableReasons(List.of(FusionReportableReason.NOT_KNOWN))
                .phased(FusionPhasedType.INFRAME)
                .likelihood(FusionLikelihoodType.NA)
                .fivePrimeVcfId("0")
                .threePrimeVcfId("1")
                .fivePrimeCoords("chr1:1")
                .threePrimeCoords("chr1:2")
                .chainLength(1)
                .chainLinks(1)
                .chainTerminated(false)
                .domainsKept("")
                .domainsLost("")
                .skippedExonsUp(0)
                .skippedExonsDown(0)
                .fusedExonUp(0)
                .fusedExonDown(0)
                .geneStart("")
                .geneContextStart("")
                .geneTranscriptStart("")
                .geneEnd("")
                .geneContextEnd("")
                .geneTranscriptEnd("")
                .junctionCopyNumber(0.0)
                .build();

        databaseAccess.writeBreakendsAndFusions(TEST_SAMPLE_ID, List.of(breakend5, breakend3), List.of(fusion));

        List<SvfusionRecord> fusionFecords = fetchTable(Tables.SVFUSION, SvfusionRecord.class);
        assertEquals(1, fusionFecords.size());

        List<SvbreakendRecord> breakendRecords = fetchTable(Tables.SVBREAKEND, SvbreakendRecord.class);
        assertEquals(2, breakendRecords.size());
    }

    private static LinxBreakend createBreakendWithId(int id)
    {
        return ImmutableLinxBreakend.builder()
                .id(id)
                .svId(0)
                .vcfId("0")
                .coords("chr1:1")
                .isStart(true)
                .gene("GENE1")
                .transcriptId("ENST0")
                .canonical(true)
                .geneOrientation("")
                .disruptive(false)
                .reportedStatus(ReportedStatus.REPORTED)
                .undisruptedCopyNumber(1)
                .regionType(TranscriptRegionType.UNKNOWN)
                .codingType(TranscriptCodingType.UNKNOWN)
                .biotype("")
                .exonUp(0)
                .exonDown(0)
                .exonicBasePhase(0)
                .nextSpliceExonRank(0)
                .nextSpliceExonPhase(0)
                .nextSpliceDistance(0)
                .totalExonCount(0)
                .build();
    }

    @Test
    public void canWriteDrivers()
    {
        LinxDriver driver = ImmutableLinxDriver.builder()
                .clusterId(0)
                .gene("GENE1")
                .eventType(DriverEventType.DEL)
                .build();

        databaseAccess.writeSvDrivers(TEST_SAMPLE_ID, List.of(driver));
        List<SvdriverRecord> records = fetchTable(Tables.SVDRIVER, SvdriverRecord.class);
        assertEquals(1, records.size());
    }

    @Test
    public void canWriteDriverCatalog()
    {
        List<DriverCatalog> drivers = List.of(
                createDriverWithType(DriverType.AMP),
                createDriverWithType(DriverType.DEL),
                createDriverWithType(DriverType.UNKNOWN)
        );

        EnumSet<DriverType> driverTypes = EnumSet.of(DriverType.AMP, DriverType.DEL);

        databaseAccess.writeLinxDriverCatalog(TEST_SAMPLE_ID, drivers, driverTypes);

        List<DrivercatalogRecord> driverCatalogRecords = fetchTable(Tables.DRIVERCATALOG, DrivercatalogRecord.class);
        assertEquals(driverTypes.size(), driverCatalogRecords.size());
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
}
