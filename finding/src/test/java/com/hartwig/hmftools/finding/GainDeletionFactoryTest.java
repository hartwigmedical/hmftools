package com.hartwig.hmftools.finding;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.finding.clinicaltranscript.ClinicalTranscriptsModelTest;
import com.hartwig.hmftools.finding.datamodel.CurationApplierTest;
import com.hartwig.hmftools.finding.datamodel.GainDeletion;
import com.hartwig.hmftools.finding.datamodel.TestFindingFactory;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFields;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFieldsBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.finding.datamodel.driver.DriverSource;
import com.hartwig.hmftools.finding.datamodel.driver.ReportedStatus;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GainDeletionFactoryTest
{

    private static final String RELEVANT_CLINICAL_GENE_COPY_NUMBERS_TSV =
            Resources.getResource("clinicalrelevantgenecopynumbers/clinical_relevant_gene_copy_numbers.tsv").getPath();
    private static final String CLINICAL_TRANSCRIPT_TSV =
            ClinicalTranscriptsModelTest.class.getClassLoader().getResource("clinicaltranscript/clinical_transcripts.tsv").getPath();

    private static DriverFields driverFields(String findingKey, ReportedStatus reportedStatus) {
        return DriverFieldsBuilder.builder().findingKey(findingKey).driverSource(DriverSource.SOMATIC).reportedStatus(reportedStatus).driverInterpretation(DriverInterpretation.HIGH).build();
    }

    @Test
    public void canExtractClicicalRelevantGeneCopyNumbers() throws IOException
    {
        List<PurpleCopyNumber> copyNumbers = createCopyNumbersList();

        GainDeletion gd = TestFindingFactory.gainDeletionBuilder()
                .gene("BRAF")
                .driver(driverFields("gd-1", ReportedStatus.CANDIDATE))
                .build();
        PurpleGeneCopyNumber clinicalRelevantGeneCopyNumbers = builder().gene("EGFR").chromosome("1").chromosomeBand("p").build();
        ArmCopyNumberFactory cn = new ArmCopyNumberFactory(
                copyNumbers, 2.0, Gender.FEMALE, OrangeRefGenomeVersion.V37);
        FindingConfig config =
                FindingConfig.createFindingConfig(Path.of(CLINICAL_TRANSCRIPT_TSV), Path.of(RELEVANT_CLINICAL_GENE_COPY_NUMBERS_TSV), null, OrangeRefGenomeVersion.V37, Gender.FEMALE, false);

        List<GainDeletion> clinicalRelevant =
                GainDeletionFactory.ClinicalRelevantGeneCopyNumber(List.of(clinicalRelevantGeneCopyNumbers), cn, config, List.of(gd));

        assertEquals(1, clinicalRelevant.size());
        assertEquals("EGFR", clinicalRelevant.get(0).gene());
        assertEquals("1", clinicalRelevant.get(0).chromosome());
        assertEquals("p", clinicalRelevant.get(0).chromosomeBand());

    }

    @Test
    public void canFilterClicicalRelevantGeneCopyNumbers2() throws IOException
    {
        List<PurpleCopyNumber> copyNumbers = createCopyNumbersList();
        GainDeletion gd = TestFindingFactory.gainDeletionBuilder()
                .gene("BRAF")
                .driver(driverFields("gd-1", ReportedStatus.CANDIDATE))
                .build();
        PurpleGeneCopyNumber clinicalRelevantGeneCopyNumbers = builder().gene("BRAF").chromosome("1").chromosomeBand("p").build();
        ArmCopyNumberFactory cn = new ArmCopyNumberFactory(
                copyNumbers, 2.0, Gender.FEMALE, OrangeRefGenomeVersion.V37);
        FindingConfig config =
                FindingConfig.createFindingConfig(Path.of(CLINICAL_TRANSCRIPT_TSV), Path.of(RELEVANT_CLINICAL_GENE_COPY_NUMBERS_TSV), null, OrangeRefGenomeVersion.V37, Gender.FEMALE, false);

        List<GainDeletion> clinicalRelevant =
                GainDeletionFactory.ClinicalRelevantGeneCopyNumber(List.of(clinicalRelevantGeneCopyNumbers), cn, config, List.of(gd));

        assertEquals(0, clinicalRelevant.size());
    }

    @NotNull
    public static List<PurpleCopyNumber> createCopyNumbersList()
    {
        List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();
        copyNumbers.add(copyNumberBuilder().chromosome("1").start(1).end(123035434).averageTumorCopyNumber(2D).build());
        copyNumbers.add(copyNumberBuilder()
                .chromosome("1")
                .start(123035435)
                .end(124035434)
                .averageTumorCopyNumber(300D)
                .build());
        copyNumbers.add(copyNumberBuilder()
                .chromosome("1")
                .start(124035435)
                .end(249250621)
                .averageTumorCopyNumber(3D)
                .build());
        return copyNumbers;
    }

    @NotNull
    public static ImmutablePurpleCopyNumber.Builder copyNumberBuilder()
    {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(Strings.EMPTY)
                .start(0)
                .end(0)
                .averageTumorCopyNumber(0D);
    }

    @NotNull
    public static ImmutablePurpleGeneCopyNumber.Builder builder()
    {
        return ImmutablePurpleGeneCopyNumber.builder()
                .gene(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .transcript(Strings.EMPTY)
                .isCanonical(true)
                .minCopyNumber(0)
                .maxCopyNumber(0)
                .minMinorAlleleCopyNumber(0);
    }
}