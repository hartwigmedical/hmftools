package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.rna.ImmutableRnaFusion;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.rna.RnaFusionFile;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.SourceData;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.TruthsetCache;
import com.hartwig.hmftools.compar.common.TruthsetValue;

import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class RnaFusionComparerTest
{
    @Rule
    public TemporaryFolder temp = new TemporaryFolder();

    private static final String SAMPLE_ID = "TEST";

    private static final RnaFusion TMPRSS2_ERG = createFusion("TMPRSS2_ERG", "21", 1000, "21", 2000, KnownFusionType.KNOWN_PAIR, 100);
    private static final RnaFusion EML4_ALK = createFusion("EML4_ALK", "2", 3000, "2", 4000, KnownFusionType.KNOWN_PAIR, 10);
    private static final RnaFusion UNKNOWN_FUSION = createFusion("GENEA_GENEB", "5", 5000, "7", 6000, KnownFusionType.NONE, 8);
    private static final RnaFusion FGFR3_TACC3 = createFusion("FGFR3_TACC3", "4", 7000, "4", 8000, KnownFusionType.KNOWN_PAIR, 12);
    private static final RnaFusion ETV6_NTRK3 = createFusion("ETV6_NTRK3", "12", 9000, "15", 10000, KnownFusionType.KNOWN_PAIR, 4);

    private static final String UNFILTERED_HEADER = String.join(TSV_DELIM,
            "KnownFusionType", "GeneNameUp", "ChrUp", "PosUp", "OrientUp", "JuncTypeUp", "CoverageUp", "MaxAnchorLengthUp",
            "GeneNameDown", "ChrDown", "PosDown", "OrientDown", "JuncTypeDown", "CoverageDown", "MaxAnchorLengthDown",
            "SVType", "SplitFrags", "RealignedFrags", "DiscordantFrags", "CohortCount", "Filter");

    private RnaFusionComparer comparer;
    private File oldDir;
    private File newDir;

    @Before
    public void setUp() throws IOException
    {
        comparer = new RnaFusionComparer(new ComparConfig(), Collections.emptyMap());
        oldDir = temp.newFolder("old");
        newDir = temp.newFolder("new");
    }

    @Test
    public void passingFileDeterminesPassStatus() throws IOException
    {
        for(SourceType source : SourceType.values())
        {
            File dir = source == SourceType.OLD ? oldDir : newDir;
            Files.write(Paths.get(RnaFusionFile.generateFilename(dir.getAbsolutePath(), SAMPLE_ID)), List.of(
                    RnaFusionFile.header() + TSV_DELIM + "Filter",
                    RnaFusionFile.write(EML4_ALK) + TSV_DELIM + "FRAGMENT_COUNT"));
            List<ComparableItem> items = comparer.loadFromFile(SAMPLE_ID, null, isofoxPaths(source, dir));

            assertEquals(1, items.size());
            RnaFusionData fusion = (RnaFusionData)items.get(0);
            assertEquals(PASS_FILTER, fusion.Filter);
            assertEquals(PASS_FILTER, fusion.fieldValues().get("Filter").displayValue());
        }
    }

    @Test
    public void oneSidedFusionsFoundInUnfilteredFile() throws IOException
    {
        writePassFusions(oldDir, List.of(TMPRSS2_ERG, EML4_ALK, UNKNOWN_FUSION));

        writeUnfilteredFusions(oldDir, Map.of(
                TMPRSS2_ERG, PASS_FILTER,
                EML4_ALK, PASS_FILTER,
                UNKNOWN_FUSION, PASS_FILTER,
                ImmutableRnaFusion.builder().from(FGFR3_TACC3).splitFragments(3).build(), "ALLELE_FREQUENCY",
                ETV6_NTRK3, "FRAGMENT_COUNT"));

        writePassFusions(newDir, List.of(TMPRSS2_ERG, FGFR3_TACC3));

        writeUnfilteredFusions(newDir, Map.of(
                TMPRSS2_ERG, PASS_FILTER,
                ImmutableRnaFusion.builder().from(EML4_ALK).splitFragments(2).build(), "FRAGMENT_COUNT",
                UNKNOWN_FUSION, PASS_FILTER,
                FGFR3_TACC3, PASS_FILTER,
                ImmutableRnaFusion.builder().from(ETV6_NTRK3).splitFragments(2).build(), "FRAGMENT_COUNT"));

        for(MatchLevel matchLevel : List.of(MatchLevel.DETAILED, MatchLevel.REPORTABLE))
        {
            List<Mismatch> mismatches = compare(
                    comparer.loadFromFile(SAMPLE_ID, null, isofoxPaths(SourceType.OLD, oldDir)),
                    comparer.loadFromFile(SAMPLE_ID, null, isofoxPaths(SourceType.NEW, newDir)), matchLevel, false);
            assertEquals(matchLevel == MatchLevel.DETAILED ? 3 : 2, mismatches.size());

            Map<String,Mismatch> mismatchesByName = mismatches.stream()
                    .collect(Collectors.toMap(x -> ((RnaFusionData)x.nonNullItem()).Name, x -> x));
            Mismatch filteredInNew = mismatchesByName.get("EML4_ALK");
            assertMismatch(filteredInNew, MismatchType.OLD_ONLY, "SplitFrags(10/2)", "Filter(PASS/FRAGMENT_COUNT)");
            assertEquals("FRAGMENT_COUNT", ((RnaFusionData)filteredInNew.NewItem).Filter);

            Mismatch filteredInOld = mismatchesByName.get("FGFR3_TACC3");
            assertMismatch(filteredInOld, MismatchType.NEW_ONLY, "SplitFrags(3/12)", "Filter(ALLELE_FREQUENCY/PASS)");
            assertEquals("ALLELE_FREQUENCY", ((RnaFusionData)filteredInOld.OldItem).Filter);

            Mismatch unknownFusion = mismatchesByName.get("GENEA_GENEB");
            if(matchLevel == MatchLevel.DETAILED)
            {
                assertMismatch(unknownFusion, MismatchType.OLD_ONLY);
                assertNull(unknownFusion.NewItem);
            }
            else
            {
                assertNull(unknownFusion);
            }

            assertNull(mismatchesByName.get("ETV6_NTRK3"));
        }
    }

    @Test
    public void truthsetFusionsArePassingInBothDirections() throws IOException
    {
        String matchedKey = "TMPRSS2:ERG:21:1000:21:2000";
        String filteredKey = "EML4:ALK:2:3000:2:4000";

        for(SourceType truthsetSource : SourceType.values())
        {
            comparer = createTruthsetComparer(truthsetSource);
            File dir = truthsetSource == SourceType.OLD ? newDir : oldDir;
            writePassFusions(dir, List.of(TMPRSS2_ERG));
            writeUnfilteredFusions(dir, Map.of(
                    TMPRSS2_ERG, PASS_FILTER,
                    ImmutableRnaFusion.builder().from(EML4_ALK).splitFragments(2).build(), "FRAGMENT_COUNT"));

            List<ComparableItem> truthsetItems = comparer.loadFromTruthset(Map.of(
                    matchedKey, List.of(new TruthsetValue(matchedKey, "SplitFrags", "100", false)),
                    filteredKey, List.of(new TruthsetValue(filteredKey, "SplitFrags", "10", false))));
            Mismatch mismatch = singleMismatch(compareTruthset(truthsetSource, truthsetItems, MatchLevel.REPORTABLE));

            if(truthsetSource == SourceType.OLD)
                assertMismatch(mismatch, MismatchType.OLD_ONLY, "SplitFrags(10/2)", "Filter(PASS/FRAGMENT_COUNT)");
            else
                assertMismatch(mismatch, MismatchType.NEW_ONLY, "SplitFrags(2/10)", "Filter(FRAGMENT_COUNT/PASS)");

            truthsetItems = comparer.loadFromTruthset(Map.of(
                    filteredKey, List.of(new TruthsetValue(filteredKey, "Filter", "FRAGMENT_COUNT", false))));
            writePassFusions(dir, List.of(EML4_ALK));

            assertEquals(0, compareTruthset(truthsetSource, truthsetItems, MatchLevel.REPORTABLE).size());
        }
    }

    @Test
    public void groupedTruthsetFusionsIncludeUnfilteredContextInBothDirections() throws IOException
    {
        for(SourceType truthsetSource : SourceType.values())
        {
            comparer = createTruthsetComparer(truthsetSource);
            String key = "EML4:ALK:2:3000:2:4000";
            List<ComparableItem> truthsetItems = comparer.loadFromTruthset(Map.of(
                    key, List.of(new TruthsetValue(key, "SplitFrags", "10", false))));

            File dir = truthsetSource == SourceType.OLD ? newDir : oldDir;
            writePassFusions(dir, List.of(ImmutableRnaFusion.builder().from(EML4_ALK).positionUp(3001).build()));
            writeUnfilteredFusions(dir, Map.of(
                    ImmutableRnaFusion.builder().from(EML4_ALK).splitFragments(2).build(), "FRAGMENT_COUNT"));

            for(boolean hasUnfilteredFile : List.of(true, false))
            {
                if(!hasUnfilteredFile)
                    Files.delete(Paths.get(RnaFusionFile.generateUnfilteredFilename(dir.getAbsolutePath(), SAMPLE_ID)));

                List<Mismatch> mismatches = compareTruthset(truthsetSource, truthsetItems, MatchLevel.DETAILED);

                String truthsetInfo = "coords(3000-4000 frags=10)";
                if(hasUnfilteredFile)
                    truthsetInfo += " unfiltered(coords(3000-4000 frags=2) filter=FRAGMENT_COUNT)";
                String fileInfo = "coords(3001-4000 frags=10)";
                assertMismatch(singleMismatch(mismatches), MismatchType.VALUE,
                        "fusions(" + (truthsetSource == SourceType.OLD ?
                                truthsetInfo + "/" + fileInfo : fileInfo + "/" + truthsetInfo) + ")");
            }
        }
    }

    @Test
    public void unavailableUnfilteredFileClearsPreviousContext() throws IOException
    {
        writePassFusions(oldDir, List.of(EML4_ALK));
        writePassFusions(newDir, Collections.emptyList());

        for(String contents : Arrays.asList(null, "", "invalid header", UNFILTERED_HEADER + "\ninvalid row"))
        {
            writeUnfilteredFusions(newDir, Map.of(EML4_ALK, "FRAGMENT_COUNT"));
            comparer.loadFromFile(SAMPLE_ID, null, isofoxPaths(SourceType.NEW, newDir));
            if(contents == null)
                Files.delete(Paths.get(RnaFusionFile.generateUnfilteredFilename(newDir.getAbsolutePath(), SAMPLE_ID)));
            else
                Files.writeString(Paths.get(RnaFusionFile.generateUnfilteredFilename(newDir.getAbsolutePath(), SAMPLE_ID)), contents);

            List<ComparableItem> oldItems = comparer.loadFromFile(SAMPLE_ID, null, isofoxPaths(SourceType.OLD, oldDir));
            List<ComparableItem> newItems = comparer.loadFromFile(SAMPLE_ID, null, isofoxPaths(SourceType.NEW, newDir));

            assertEquals(0, newItems.size());

            List<Mismatch> mismatches = compare(oldItems, newItems, MatchLevel.DETAILED, false);

            Mismatch mismatch = singleMismatch(mismatches);
            assertMismatch(mismatch, MismatchType.OLD_ONLY);
            assertNull(mismatch.NewItem);
        }
    }

    @Test
    public void unfilteredFusionCoordinatesAreLifted() throws IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        ComparConfig.addConfig(configBuilder);
        configBuilder.setValue("sample", SAMPLE_ID);
        configBuilder.setValue("sample_dir_old", oldDir.getAbsolutePath());
        configBuilder.setValue("sample_dir_new", newDir.getAbsolutePath());
        configBuilder.setValue(ComparConfig.REQUIRES_LIFTOVER, true);
        ComparConfig config = new ComparConfig(configBuilder);
        config.LiftoverCache.getChromosomeMappings("2", true).clear();
        config.LiftoverCache.addMapping("2", 1, 10000, 101, 10100, false);
        comparer = new RnaFusionComparer(config, Collections.emptyMap());

        writePassFusions(oldDir, Collections.emptyList());
        writeUnfilteredFusions(oldDir, Map.of(EML4_ALK, "FRAGMENT_COUNT"));
        RnaFusion liftedFusion = ImmutableRnaFusion.builder().from(EML4_ALK)
                .chromosomeUp("chr2").chromosomeDown("chr2").positionUp(3100).positionDown(4100).build();
        writePassFusions(newDir, List.of(liftedFusion));

        List<Mismatch> mismatches = compare(
                comparer.loadFromFile(SAMPLE_ID, null, isofoxPaths(SourceType.OLD, oldDir)),
                comparer.loadFromFile(SAMPLE_ID, null, isofoxPaths(SourceType.NEW, newDir)), MatchLevel.DETAILED, false);

        Mismatch mismatch = singleMismatch(mismatches);
        assertMismatch(mismatch, MismatchType.NEW_ONLY, "Filter(FRAGMENT_COUNT/PASS)");
        RnaFusionData counterpart = (RnaFusionData)mismatch.OldItem;
        assertEquals(liftedFusion.positionUp(), counterpart.PositionUp);
        assertEquals(liftedFusion.positionDown(), counterpart.PositionDown);
    }

    private static RnaFusionComparer createTruthsetComparer(final SourceType truthsetSource)
    {
        ComparConfig config = new ComparConfig();
        config.Sources.set(truthsetSource == SourceType.OLD ? 0 : 1,
                SourceData.fromTruthsetSource(truthsetSource, new TruthsetCache()));
        return new RnaFusionComparer(config, Collections.emptyMap());
    }

    private List<Mismatch> compareTruthset(
            final SourceType truthsetSource, final List<ComparableItem> truthsetItems, final MatchLevel matchLevel)
    {
        boolean truthsetIsOld = truthsetSource == SourceType.OLD;
        List<ComparableItem> fileItems = comparer.loadFromFile(SAMPLE_ID, null,
                truthsetIsOld ? isofoxPaths(SourceType.NEW, newDir) : isofoxPaths(SourceType.OLD, oldDir));
        return compare(truthsetIsOld ? truthsetItems : fileItems, truthsetIsOld ? fileItems : truthsetItems, matchLevel, true);
    }

    private List<Mismatch> compare(
            final List<ComparableItem> oldItems, final List<ComparableItem> newItems, final MatchLevel matchLevel,
            boolean includesTruthset)
    {
        List<Mismatch> mismatches = Lists.newArrayList();
        comparer.compareItems(mismatches, matchLevel, false, includesTruthset, oldItems, newItems);
        return mismatches;
    }

    private static Mismatch singleMismatch(final List<Mismatch> mismatches)
    {
        assertEquals(1, mismatches.size());
        return mismatches.get(0);
    }

    private static void assertMismatch(final Mismatch mismatch, final MismatchType type, final String... differences)
    {
        assertNotNull(mismatch);
        assertEquals(type, mismatch.Type);
        assertEquals(Arrays.asList(differences), mismatch.DiffValues);
    }

    private static RnaFusion createFusion(
            final String name, final String chromosomeUp, final int positionUp, final String chromosomeDown, final int positionDown,
            final KnownFusionType knownType, final int splitFragments)
    {
        return ImmutableRnaFusion.builder()
                .name(name)
                .knownType(knownType)
                .chromosomeUp(chromosomeUp)
                .chromosomeDown(chromosomeDown)
                .positionUp(positionUp)
                .positionDown(positionDown)
                .orientationUp((byte)1)
                .orientationDown((byte)-1)
                .junctionTypeUp("KNOWN")
                .junctionTypeDown("CANONICAL")
                .transcriptUp("ENST01")
                .transcriptDown("ENST02")
                .exonUp(2)
                .exonDown(4)
                .svType(StructuralVariantType.BND)
                .splitFragments(splitFragments)
                .realignedFrags(1)
                .discordantFrags(2)
                .depthUp(200)
                .depthDown(300)
                .maxAnchorLengthUp(40)
                .maxAnchorLengthDown(50)
                .cohortFrequency(0)
                .build();
    }

    private static PipelineSourcePaths isofoxPaths(final SourceType sourceType, final File isofoxDir)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        PipelineSourcePaths.registerConfig(configBuilder);
        configBuilder.setValue(ISOFOX_DIR_CFG + "_" + sourceType.configStr(), isofoxDir.getAbsolutePath());
        return PipelineSourcePaths.fromConfig(sourceType, configBuilder);
    }

    private static void writePassFusions(final File isofoxDir, final List<RnaFusion> fusions) throws IOException
    {
        List<String> lines = Lists.newArrayList(RnaFusionFile.header());
        fusions.forEach(x -> lines.add(RnaFusionFile.write(x)));
        Files.write(Paths.get(RnaFusionFile.generateFilename(isofoxDir.getAbsolutePath(), SAMPLE_ID)), lines);
    }

    private static void writeUnfilteredFusions(final File isofoxDir, final Map<RnaFusion,String> fusionFilters) throws IOException
    {
        List<String> lines = Lists.newArrayList(UNFILTERED_HEADER);

        for(Map.Entry<RnaFusion,String> entry : fusionFilters.entrySet())
        {
            RnaFusion fusion = entry.getKey();

            lines.add(String.join(TSV_DELIM,
                    fusion.knownType().toString(), RnaFusionFile.geneUp(fusion), fusion.chromosomeUp(),
                    String.valueOf(fusion.positionUp()), String.valueOf(fusion.orientationUp()), fusion.junctionTypeUp(),
                    String.valueOf(fusion.depthUp()), String.valueOf(fusion.maxAnchorLengthUp()),
                    RnaFusionFile.geneDown(fusion), fusion.chromosomeDown(), String.valueOf(fusion.positionDown()),
                    String.valueOf(fusion.orientationDown()), fusion.junctionTypeDown(), String.valueOf(fusion.depthDown()),
                    String.valueOf(fusion.maxAnchorLengthDown()), fusion.svType().toString(), String.valueOf(fusion.splitFragments()),
                    String.valueOf(fusion.realignedFrags()), String.valueOf(fusion.discordantFrags()),
                    String.valueOf(fusion.cohortFrequency()), entry.getValue()));
        }

        Files.write(Paths.get(RnaFusionFile.generateUnfilteredFilename(isofoxDir.getAbsolutePath(), SAMPLE_ID)), lines);
    }
}
