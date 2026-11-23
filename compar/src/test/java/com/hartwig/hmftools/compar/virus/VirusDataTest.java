package com.hartwig.hmftools.compar.virus;

import static com.hartwig.hmftools.compar.ComparTestUtil.assertDifferencesAreForFields;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.viridian.ClonalityStatus;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.common.virusbreakend.VirusBreakendQCStatus;
import com.hartwig.hmftools.common.virusinterpreter.ImmutableAnnotatedVirus;
import com.hartwig.hmftools.common.virusinterpreter.VirusLikelihoodType;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ComparableItemTest;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;
import com.hartwig.hmftools.compar.virus.VirusComparer.Fields;

import org.jetbrains.annotations.Nullable;
import org.junit.Before;
import org.junit.Test;

public class VirusDataTest extends ComparableItemTest<VirusData, VirusComparer, TestVirusDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new VirusComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestVirusDataBuilder.BUILDER;

        fieldToAlternateValueInitializer = Map.of(
                Fields.VirusName.toString(), b -> b.virusName = "Human papillomavirus type 18",
                Fields.Interpretation.toString(), b -> b.reportingType = VirusType.EBV,
                Fields.Integrations.toString(), b -> b.integrations = 40,
                Fields.MeanCoverage.toString(), b -> b.meanDepth = 2,
                Fields.PercentageCovered.toString(), b -> b.coverage = 0.2,
                Fields.CopiesPerTumorCell.toString(), b -> b.copiesPerTumorCell = 0.2,
                Fields.ClonalityStatus.toString(), b -> b.clonalityStatus = ClonalityStatus.NOT_CLONAL);
        nameToAlternateIndexInitializer = Map.of("virusGroup", b -> b.virusGroup = "Human gammaherpesvirus 4");
        reportabilityFieldToFalseReportabilityInitializer = Map.of(Fields.Reported.toString(), b -> b.reported = false);
        nameToNonPassInitializer = Collections.emptyMap();
    }

    @Test
    public void testFindMismatchVirusInterpreterToViridian()
    {
        // VirusInterpreter names the reference genome's taxon, writes coverage as a percentage, and has no copies per tumor cell
        // or clonality status, so only the name differs
        VirusData virusInterpreterData = virusInterpreterData("Human papillomavirus 16", VirusType.HPV);
        VirusData viridianData = builder.create();

        assertTrue(virusInterpreterData.matches(viridianData));
        Mismatch mismatch = virusInterpreterData.findMismatch(comparer, viridianData, MatchLevel.DETAILED, false);
        assertDifferencesAreForFields(Set.of(Fields.VirusName.toString()), mismatch.DiffValues);
    }

    @Test
    public void testFromVirusInterpreterKey()
    {
        List<VirusData> data = List.of(
                virusInterpreterData("Human papillomavirus 16", VirusType.HPV),
                virusInterpreterData("human gammaherpesvirus 4", VirusType.EBV),
                virusInterpreterData("Human gammaherpesvirus 8", VirusType.HHV8),
                virusInterpreterData("hepatitis B virus genotype A", VirusType.HBV),
                virusInterpreterData("Merkel cell polyomavirus", VirusType.MCV),
                virusInterpreterData("Hepatitis C virus subtype 1b", null),
                virusInterpreterData("hepatitis C virus genotype 1a", null),
                virusInterpreterData("Hepacivirus hominis", null),
                virusInterpreterData("Betapolyomavirus hominis", null),
                virusInterpreterData("Human T-cell leukemia virus type I", null),
                virusInterpreterData("JC polyomavirus", null),
                virusInterpreterData("Human papillomavirus type 5", null));

        List<String> expectedKeys = List.of(
                "Human papillomavirus",
                "Human gammaherpesvirus 4",
                "Human gammaherpesvirus 8",
                "Hepatitis B virus",
                "Merkel cell polyomavirus",
                "Hepatitis C virus",
                "Hepatitis C virus",
                "Hepatitis C virus",
                "BK polyomavirus",
                "Human T-lymphotropic virus 1",
                "JC polyomavirus",
                "Human papillomavirus type 5");

        assertEquals(expectedKeys, data.stream().map(VirusData::key).toList());
    }

    @Test
    public void testDuplicateKeys()
    {
        // two HPV types in one VirusInterpreter output map to the same virus group
        List<ComparableItem> items = List.of(
                virusInterpreterData("Human papillomavirus 16", VirusType.HPV),
                virusInterpreterData("human papillomavirus 18", VirusType.HPV),
                virusInterpreterData("human gammaherpesvirus 4", VirusType.EBV));

        assertEquals(List.of("Human papillomavirus"), VirusComparer.duplicateKeys(items));
    }

    @Test
    public void testCompareItemsDuplicateKeys()
    {
        // the first VirusInterpreter row in the virus group pairs with Viridian's row and the second shows as old only
        VirusData hpv16 = virusInterpreterData("Human papillomavirus 16", VirusType.HPV);
        VirusData hpv18 = virusInterpreterData("human papillomavirus 18", VirusType.HPV);
        VirusData viridianData = builder.create();

        List<Mismatch> mismatches = new ArrayList<>();
        comparer.compareItems(
                mismatches, MatchLevel.DETAILED, false, false, new ArrayList<>(List.of(hpv16, hpv18)),
                new ArrayList<>(List.of(viridianData)));

        List<Mismatch> expected = List.of(
                new Mismatch(hpv16, viridianData, MismatchType.VALUE, Collections.emptyList()),
                new Mismatch(hpv18, null, MismatchType.OLD_ONLY, Collections.emptyList()));
        assertEquals(expected, mismatches);
    }

    @Test
    public void testFindMismatchPurpleFailure()
    {
        // without a usable Purple fit Viridian writes no copies per tumor cell, so only the clonality status differs
        VirusData purpleFailData = builder.create(b ->
        {
            b.copiesPerTumorCell = null;
            b.clonalityStatus = ClonalityStatus.PURPLE_FAIL_NO_TUMOR;
        });
        VirusData purplePassData = builder.create();

        Mismatch mismatch = purplePassData.findMismatch(comparer, purpleFailData, MatchLevel.DETAILED, false);
        assertDifferencesAreForFields(Set.of(Fields.ClonalityStatus.toString()), mismatch.DiffValues);
    }

    private VirusData virusInterpreterData(final String name, @Nullable final VirusType interpretation)
    {
        return VirusData.fromVirusInterpreter(
                ImmutableAnnotatedVirus.builder()
                        .taxid(1)
                        .name(name)
                        .qcStatus(VirusBreakendQCStatus.NO_ABNORMALITIES)
                        .integrations(4)
                        .interpretation(interpretation)
                        .percentageCovered(90)
                        .meanCoverage(20)
                        .expectedClonalCoverage(5.0)
                        .reported(true)
                        .blacklisted(false)
                        .virusDriverLikelihoodType(VirusLikelihoodType.HIGH)
                        .build(),
                comparer.fieldsList());
    }
}
